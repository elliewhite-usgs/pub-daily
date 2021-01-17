# measures of fit with a modified bR2
goffuncs <- list.files("libraries/hydroGOFm/R")
for(i in 1:length(goffuncs)){
  source(paste0("Libraries/hydroGOFm/R/", goffuncs[i]))
}
remove(goffuncs)

# output of the make array function is in dimensions of [nbasins, tsteps, nfeat]
makearray <- function(tset, tsetpvs){ 
  tset_long <- melt(tset, id.vars=c("CDEC_ID", "DATE"), measure.vars=c(colnames(tsetpvs), "FLOW"), variable.name="VARS", value.name="VALUE")
  tset_list <- split(tset_long, tset_long$VARS)
  tset_list <- lapply(tset_list, function(x) x[!(names(x) %in% "VARS")]) # drop the VARS column
  tset_list <- lapply(tset_list, function(x) dcast(x, CDEC_ID~DATE, value.var = "VALUE")) # make dates as columns
  tset_list <- lapply(tset_list, function(x) x[!(names(x) %in% "CDEC_ID")]) # drop the CDEC_ID column
  tset_array <- array(unlist(tset_list), dim=c(dim(tset_list[[1]]), length(tset_list)))
}

# the 3D output of the makearray function is turned into a sliding window for the LSTM in dimensions of [nbasins*(tsteps-ncells+1), ncells, nfeat]
reshapearray <- function(tset_array, ncells){
  nbasins <- dim(tset_array)[1]
  tsteps <- dim(tset_array)[2]
  nfeat <- dim(tset_array)[3]
  new_tset_array <- array(NA, dim=c(nbasins*(tsteps-ncells+1), ncells, nfeat))
  for(i in 1:nbasins){
    for (j in 1:(tsteps-ncells+1)){
      for (n in 1:ncells){
        for (k in 1:nfeat){
          new_tset_array[(i-1)*(tsteps-ncells+1)+j, n, k] <- tset_array[i, j+n-1, k]
        }      
      }
    }
  }
  return(new_tset_array)
}

# trains the model with various parameters
# grouping style can be "logo" or "lmgo"
# activation_h is the hidden layer activation, make it nonlinear
# activation_l is the last layer activation, make it match your response variable
# dropout rate is set to 0 by default. If set it will drop out nodes after the LSTM layer and reduce by half and dropout nodes after the hidden layer 
# define a factor to do reduce_lr_on_plateau
# define a patience_early to do early stopping
# ncells is the size of the moving window or connected sequence of observations for LSTM
# set verbose to 0 when ready to run fast 

lstmcv <- function(data, groupingstyle, lstm_units=64, first_units=64, hidden_units=64, output_units=1, activation_h="tanh", activation_l="exponential", nepochs=100, dropout_rate=0, ncells=5, factor=NULL, patience_early=NULL, verbose=1){
  nbasins <- length(unique(data$CDEC_ID))
  tsteps <- aggregate(FLOW~CDEC_ID, data=data, FUN=length)[1, 2] # should be the same for all basins
  nfeat <- dim(data[,c(10,13:(ncol(data)-1))])[2]

  nnmodel <- keras_model_sequential() %>%
    layer_lstm(units = lstm_units, input_shape = c(ncells, nfeat), return_sequences = FALSE) %>% 
    layer_dropout(rate = dropout_rate) %>% 
    layer_dense(units = first_units, activation = activation_h) %>% 
    layer_dropout(rate = dropout_rate/2) %>%
    layer_dense(units = hidden_units, activation = activation_l) %>% 
    layer_dense(units = output_units)
  
  nnmodel %>% 
    compile(optimizer = optimizer_adam(clipnorm = TRUE, clipvalue = 1), loss = loss_mean_squared_error, metrics = c("mae"))
  
  if(!is.null(factor)){
    callback1 <- callback_reduce_lr_on_plateau(monitor = "val_loss", patience=100, factor=factor)
    if(!is.null(patience_early)){
      callback2 <- callback_early_stopping(mode="min", restore_best_weights=TRUE, patience=patience_early)
      callbacks <- list(callback1, callback2)
    } else {
      callbacks <- callback1
    }
  } else {
    if(!is.null(patience_early)){
      callbacks <- callback_early_stopping(mode="min", restore_best_weights=TRUE, patience=patience_early)
    } else {
      callbacks <- NULL
    }
  }
  
  # for leave one group (basin) out cross validation
  if(groupingstyle=="logo"){
    results <- modgof <- list()
    
    for (k in 1:(length(unique(data$CDEC_ID)))){
      h <- unique(data$CDEC_ID)[k]
      trainset <- data[data$CDEC_ID!=h,]
      testset <- data[data$CDEC_ID==h,]
      
      trainsetpvs <- trainset[,c(10,13:(ncol(trainset)-1))]
      testsetpvs <- testset[,c(10,13:(ncol(testset)-1))]
      
      # mean-standard deviation normalization
      ## find mean and sd column-wise of training data
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      
      ## to just center: sweep(trainsetpvs, 2L, trainmean, FUN="-")
      ## centered and scaled, note that testset is centered by trainset mean and sd (no leakage)!
      trainsetpvs_normalized <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs_normalized <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      # make new dataframes with the normalized predictor variables
      trainset <- cbind(CDEC_ID=trainset$CDEC_ID, DATE=trainset$DATE, trainsetpvs_normalized, FLOW=trainset$FLOW)
      testset <- cbind(CDEC_ID=testset$CDEC_ID, DATE=testset$DATE, testsetpvs_normalized, FLOW=testset$FLOW)
      
      # make matrices, note that only predictor variables get normalized
      trainsetpvs <- as.matrix(trainsetpvs_normalized)
      testsetpvs <- as.matrix(testsetpvs_normalized)
      trainsetrv <- as.matrix(trainset$FLOW)
      testsetrv <- as.matrix(testset$FLOW)
      
      # make an array out of matrices
      trainset_array <- makearray(trainset, trainsetpvs)
      testset_array <- makearray(testset, testsetpvs)
      
      # reshape arrays into moving window for LSTM training
      trainset_array_window <- reshapearray(trainset_array, ncells)
      testset_array_window <- reshapearray(testset_array, ncells)
      
      # split into predictor variables and response variables
      trainsetpv_array <- trainset_array_window[ , , 1:nfeat]
      trainsetrv_array <- trainset_array_window[ , , nfeat+1]
      testsetpv_array <- testset_array_window[ , , 1:nfeat]
      testsetrv_array <- testset_array_window[ , , nfeat+1]
      
      # reset lr in case a callback in a previous loop has changed it
      lr <- 1e-3
      
      nnmodel %>%
        fit(trainsetpv_array, trainsetrv_array, epochs = nepochs, verbose=verbose, batch_size=320, validation_split = 0.2, shuffle = FALSE, callbacks = callbacks) 
      
      predictions <- nnmodel %>% predict(testsetpv_array)
      predictions <- predictions[ , 1]
      results[[k]] <- cbind(obs=testsetrv_array[, ncells], pred=predictions)
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
    }
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
    colnames(modgof) <- unique(data$CDEC_ID)
  }
  
  # for random 5-fold leave multiple groups out cross validation, meaning a random 1/5 of the basins will be left out of training
  if(groupingstyle=="lmgo"){
    nfolds <- 5
    group <- as.data.frame(unique(data$CDEC_ID))
    colnames(group) <- "CDEC_ID"
    group$kftrain <- kfold(nrow(group), nfolds)
    data <- merge(data, group, by="CDEC_ID")
    
    results <- modgof <- list()
    for(k in 1:nfolds) {
      trainset <- data[data$kftrain != k, ]
      testset <- data[data$kftrain == k, ]
      
      trainsetpvs <- trainset[,c(10,13:(ncol(trainset)-1))]
      testsetpvs <- testset[,c(10,13:(ncol(testset)-1))]
      
      # mean-standard deviation normalization
      ## find mean and sd column-wise of training data
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      
      ## to just center: sweep(trainsetpvs, 2L, trainmean, FUN="-")
      ## centered and scaled, note that testset is centered by trainset mean and sd (no leakage)!
      trainsetpvs_normalized <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs_normalized <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      # make new dataframes with the normalized predictor variables
      trainset <- cbind(CDEC_ID=trainset$CDEC_ID, DATE=trainset$DATE, trainsetpvs_normalized, FLOW=trainset$FLOW)
      testset <- cbind(CDEC_ID=testset$CDEC_ID, DATE=testset$DATE, testsetpvs_normalized, FLOW=testset$FLOW)
      
      # make matrices, note that only predictor variables get normalized
      trainsetpvs <- as.matrix(trainsetpvs_normalized)
      testsetpvs <- as.matrix(testsetpvs_normalized)
      trainsetrv <- as.matrix(trainset$FLOW)
      testsetrv <- as.matrix(testset$FLOW)
      
      # make an array out of matrices
      trainset_array <- makearray(trainset, trainsetpvs)
      testset_array <- makearray(testset, testsetpvs)
      
      # reshape arrays into moving window for LSTM training
      trainset_array_window <- reshapearray(trainset_array, ncells)
      testset_array_window <- reshapearray(testset_array, ncells)
      
      # split into predictor variables and response variables
      trainsetpv_array <- trainset_array_window[ , , 1:nfeat]
      trainsetrv_array <- trainset_array_window[ , , nfeat+1]
      testsetpv_array <- testset_array_window[ , , 1:nfeat]
      testsetrv_array <- testset_array_window[ , , nfeat+1]
      
      # reset lr in case a callback in a previous loop has changed it
      lr <- 1e-3
      
      nnmodel %>%
        fit(trainsetpv_array, trainsetrv_array, epochs = nepochs, verbose=verbose, batch_size=320, validation_split = 0.2, shuffle = FALSE, callbacks = callbacks)
      
      predictions <- nnmodel %>% predict(testsetpv_array)
      predictions <- predictions[ , 1]
      results[[k]] <- cbind(obs=testsetrv_array[, ncells], pred=predictions)
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
    }
    
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
  }
  
  if(groupingstyle=="logo"){
    list(results=results, gof=modgof)
  } else if(groupingstyle=="lmgo"){
    list(results=results, gof=modgof, kf=group)
  } else{
    print("not a valid groupingstyle!")
  }
}

nnscv <- function(data, groupingstyle){
  nfeat <- dim(data[,c(10,13:(ncol(data)-1))])[2]
  
  # define the model
  nnmodel <- keras_model_sequential() %>%
    layer_dense(units=64, activation="tanh", input_shape=nfeat) %>%
    layer_dense(units=64, activation="exponential") %>%
    layer_dense(units=1)
  
  # can define with relu activations
  # nnmodel <- keras_model_sequential() %>%
  #   layer_dense(units=64, activation="relu", input_shape=nfeat) %>%
  #   layer_dense(units=64, activation="relu") %>%
  #   layer_dense(units=1)
  
  # # can define with dropout, but it doesn't help the error, maybe because it needs larger units
  # nnmodel <- keras_model_sequential() %>%
  #   layer_dense(units=64+13, activation="tanh", input_shape=nfeat) %>%
  #   layer_dropout(rate = 0.20) %>% 
  #   layer_dense(units=64+7, activation="exponential") %>%
  #   layer_dropout(rate = 0.10) %>% 
  #   layer_dense(units=1)
  
  # compile 
  nnmodel %>%
    compile(optimizer=optimizer_adam(clipnorm = TRUE, clipvalue = 1), loss=loss_mean_squared_error, metrics=c("mae"))
  
  # define a callback for early stopping
  callbacks <- list(callback_early_stopping(patience=10, mode="min", restore_best_weights=TRUE)) # if you add this you have to make sure lr resets after every loop!!! callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.1, patience=10)
  
  # for resubstitution, where testing data and training data is the same-----------------------------------------
  if(groupingstyle=="resub"){
    trainsetpvs <- as.matrix(data[,c(10,13:(ncol(data)-1))])
    trainsetrv <- as.matrix(data$FLOW)
    
    # mean-standard deviation normalization. First, find mean and sd column-wise of training data
    trainmean <- apply(trainsetpvs, 2, mean)
    trainsd <- apply(trainsetpvs, 2, sd)
    
    # to just center: sweep(trainsetpvs, 2L, trainmean, FUN="-"), but we want to center and normalize by sd
    trainsetpvs <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
    
    nnmodel %>% 
      fit(trainsetpvs, trainsetrv, epochs=30, batch_size=32, verbose=1, validation_split=0.2) #callbacks=callbacks
    
    # make predictions 
    predictions <- nnmodel %>% predict(trainsetpvs)
    predictions <- predictions[ , 1] # because output layer was specified to be of unit=1
    results <- cbind(obs=trainsetrv[ , 1], pred=predictions)
    modgof <- gof(as.data.frame(results)$pred, as.data.frame(results)$obs)
    print(paste0("resub modgof=", modgof))
  }
  
  # for k-fold cross validation---------------------------------------------------------------------------------
  if(is.integer(groupingstyle)){
    group <- kfold(data, groupingstyle)
    results <- modgof <- list()
    for(k in 1:groupingstyle) {
      trainset <- data[group != k, ]
      testset <- data[group == k, ]
      
      # separate into predictor variables and response variable
      testsetpvs <- as.matrix(testset[, c(10,13:(ncol(testset)-1))])
      trainsetpvs <- as.matrix(trainset[, c(10,13:(ncol(trainset)-1))])
      testsetrv <- as.matrix(testset$FLOW)
      trainsetrv <- as.matrix(trainset$FLOW)
      
      # mean-standard deviation normalization. Note that testset is centered by trainset mean and sd (no leakage)!
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      trainsetpvs <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      nnmodel %>%
        fit(trainsetpvs, trainsetrv, epochs=30, batch_size=32, verbose=1, validation_split=0.2)
      
      predictions <- nnmodel %>% predict(testsetpvs)
      predictions <- predictions[ , 1] 
      results[[k]] <- cbind(obs=testsetrv[ , 1], pred=predictions, nforstitch=as.numeric(rownames(testset))) # adding rownames here for later stitching the dataset back to its original shape
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
      print(paste0("kfold k=", k, ", modgof=", modgof[[k]]["bR2",]))
    }
    
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
  }
  
  # for leave one group (basin) out cross validation------------------------------------------------------------
  if(groupingstyle=="logo"){
    results <- modgof <- list()
    
    for (k in 1:(length(unique(data$CDEC_ID)))){
      h <- unique(data$CDEC_ID)[k]
      testset <- data[data$CDEC_ID==h,]
      trainset <- data[data$CDEC_ID!=h,]
      testsetpvs <- as.matrix(testset[,c(10,13:(ncol(testset)-1))])
      trainsetpvs <- as.matrix(trainset[,c(10,13:(ncol(trainset)-1))])
      testsetrv <- as.matrix(testset$FLOW)
      trainsetrv <- as.matrix(trainset$FLOW)
      
      # mean-standard deviation normalization. Note that testset is centered by trainset mean and sd (no leakage)!
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      trainsetpvs <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      nnmodel %>%
        fit(trainsetpvs, trainsetrv, epochs=100, batch_size=32, verbose=1, validation_split=0.2, callbacks=callbacks)
      
      predictions <- nnmodel %>% predict(testsetpvs)
      predictions <- predictions[ , 1] # because output layer was specified to be of unit=1
      results[[k]] <- cbind(obs=testsetrv[ , 1], pred=predictions)
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
      print(paste0("logo k=", k, ", h=", h, ", modgof=", modgof[[k]]["bR2",]))
    }
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
    colnames(modgof) <- unique(data$CDEC_ID)
  }
  
  # for random 5-fold leave multiple groups out cross validation, meaning a random 1/5 of the basins will be left out of training-------------------------------------------------------------------------------------------------
  if(groupingstyle=="lmgo"){
    nfolds <- 5
    group <- as.data.frame(unique(data$CDEC_ID))
    colnames(group) <- "CDEC_ID"
    group$kftrain <- kfold(nrow(group), nfolds)
    data <- merge(data, group, by="CDEC_ID")
    
    results <- modgof <- list()
    for(k in 1:nfolds) {
      testset <- data[data$kftrain == k, ]
      trainset <- data[data$kftrain != k, ]
      testsetpvs <- as.matrix(testset[,c(10,13:(ncol(testset)-2))])
      trainsetpvs <- as.matrix(trainset[,c(10,13:(ncol(trainset)-2))])
      testsetrv <- as.matrix(testset$FLOW)
      trainsetrv <- as.matrix(trainset$FLOW)
      
      # mean-standard deviation normalization. Note that testset is centered by trainset mean and sd (no leakage)!
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      trainsetpvs <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      nnmodel %>%
        fit(trainsetpvs, trainsetrv, epochs=30, batch_size=32, verbose=1, validation_split=0.2)
      
      predictions <- nnmodel %>% predict(testsetpvs)
      predictions <- predictions[ , 1] # because output layer was specified to be of unit=1
      results[[k]] <- cbind(obs=testsetrv[ , 1], pred=predictions, nforstitch=as.numeric(rownames(testset)))
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
      print(paste0("lmgo k=", k, ", modgof=", modgof[[k]]["bR2",]))
    }
    
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
  }
  
  if(groupingstyle=="resub"){
    list(results=results, gof=modgof)
  } else if(groupingstyle=="logo"){
    list(results=results, gof=modgof)
  } else{
    list(results=results, gof=modgof, kf=group)
  }
}

# trying without the 2 layers after the lstm
lstmpcv <- function(data, groupingstyle, lstm_units=64, output_units=1, activation_h="tanh", activation_l="exponential", nepochs=100, dropout_rate=0, ncells=5, factor=NULL, patience_early=NULL, verbose=1){
  nbasins <- length(unique(data$CDEC_ID))
  tsteps <- aggregate(FLOW~CDEC_ID, data=data, FUN=length)[1, 2] # should be the same for all basins
  nfeat <- dim(data[,c(10,13:(ncol(data)-1))])[2]
  
  nnmodel <- keras_model_sequential() %>%
    layer_lstm(units = lstm_units, input_shape = c(ncells, nfeat), return_sequences = FALSE) %>% 
    layer_dropout(rate = dropout_rate) %>% 
    layer_dense(units = output_units)
  
  nnmodel %>% 
    compile(optimizer = optimizer_adam(clipnorm = TRUE, clipvalue = 1), loss = loss_mean_squared_error, metrics = c("mae"))
  
  # # add another callback so you can vizualize runs 
  # callback_tensorboard(log_dir = "logs/run_b")
  
  if(!is.null(factor)){
    callback1 <- callback_reduce_lr_on_plateau(monitor = "val_loss", patience=100, factor=factor)
    if(!is.null(patience_early)){
      callback2 <- callback_early_stopping(mode="min", restore_best_weights=TRUE, patience=patience_early)
      callbacks <- list(callback1, callback2)
    } else {
      callbacks <- callback1
    }
  } else {
    if(!is.null(patience_early)){
      callbacks <- callback_early_stopping(mode="min", restore_best_weights=TRUE, patience=patience_early)
    } else {
      callbacks <- NULL
    }
  }
  
  # for leave one group (basin) out cross validation
  if(groupingstyle=="logo"){
    results <- modgof <- list()
    
    for (k in 1:(length(unique(data$CDEC_ID)))){
      h <- unique(data$CDEC_ID)[k]
      trainset <- data[data$CDEC_ID!=h,]
      testset <- data[data$CDEC_ID==h,]
      
      trainsetpvs <- trainset[,c(10,13:(ncol(trainset)-1))]
      testsetpvs <- testset[,c(10,13:(ncol(testset)-1))]
      
      # mean-standard deviation normalization
      ## find mean and sd column-wise of training data
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      
      ## to just center: sweep(trainsetpvs, 2L, trainmean, FUN="-")
      ## centered and scaled, note that testset is centered by trainset mean and sd (no leakage)!
      trainsetpvs_normalized <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs_normalized <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      # make new dataframes with the normalized predictor variables
      trainset <- cbind(CDEC_ID=trainset$CDEC_ID, DATE=trainset$DATE, trainsetpvs_normalized, FLOW=trainset$FLOW)
      testset <- cbind(CDEC_ID=testset$CDEC_ID, DATE=testset$DATE, testsetpvs_normalized, FLOW=testset$FLOW)
      
      # make matrices, note that only predictor variables get normalized
      trainsetpvs <- as.matrix(trainsetpvs_normalized)
      testsetpvs <- as.matrix(testsetpvs_normalized)
      trainsetrv <- as.matrix(trainset$FLOW)
      testsetrv <- as.matrix(testset$FLOW)
      
      # make an array out of matrices
      trainset_array <- makearray(trainset, trainsetpvs)
      testset_array <- makearray(testset, testsetpvs)
      
      # reshape arrays into moving window for LSTM training
      trainset_array_window <- reshapearray(trainset_array, ncells)
      testset_array_window <- reshapearray(testset_array, ncells)
      
      # split into predictor variables and response variables
      trainsetpv_array <- trainset_array_window[ , , 1:nfeat]
      trainsetrv_array <- trainset_array_window[ , , nfeat+1]
      testsetpv_array <- testset_array_window[ , , 1:nfeat]
      testsetrv_array <- testset_array_window[ , , nfeat+1]
      
      # reset lr in case a callback in a previous loop has changed it
      lr <- 1e-3
      
      nnmodel %>%
        fit(trainsetpv_array, trainsetrv_array, epochs = nepochs, verbose=verbose, batch_size=320, validation_split = 0.2, shuffle = FALSE, callbacks = callbacks) 
      
      predictions <- nnmodel %>% predict(testsetpv_array)
      predictions <- predictions[ , 1]
      results[[k]] <- cbind(obs=testsetrv_array[, ncells], pred=predictions)
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
    }
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
    colnames(modgof) <- unique(data$CDEC_ID)
  }
  
  # for random 5-fold leave multiple groups out cross validation, meaning a random 1/5 of the basins will be left out of training
  if(groupingstyle=="lmgo"){
    nfolds <- 5
    group <- as.data.frame(unique(data$CDEC_ID))
    colnames(group) <- "CDEC_ID"
    group$kftrain <- kfold(nrow(group), nfolds)
    data <- merge(data, group, by="CDEC_ID")
    
    results <- modgof <- list()
    for(k in 1:nfolds) {
      trainset <- data[data$kftrain != k, ]
      testset <- data[data$kftrain == k, ]
      
      trainsetpvs <- trainset[,c(10,13:(ncol(trainset)-1))]
      testsetpvs <- testset[,c(10,13:(ncol(testset)-1))]
      
      # mean-standard deviation normalization
      ## find mean and sd column-wise of training data
      trainmean <- apply(trainsetpvs, 2, mean)
      trainsd <- apply(trainsetpvs, 2, sd)
      
      ## to just center: sweep(trainsetpvs, 2L, trainmean, FUN="-")
      ## centered and scaled, note that testset is centered by trainset mean and sd (no leakage)!
      trainsetpvs_normalized <- sweep(sweep(trainsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/") 
      testsetpvs_normalized <- sweep(sweep(testsetpvs, 2L, trainmean, FUN="-"), 2, trainsd, FUN="/")
      
      # make new dataframes with the normalized predictor variables
      trainset <- cbind(CDEC_ID=trainset$CDEC_ID, DATE=trainset$DATE, trainsetpvs_normalized, FLOW=trainset$FLOW)
      testset <- cbind(CDEC_ID=testset$CDEC_ID, DATE=testset$DATE, testsetpvs_normalized, FLOW=testset$FLOW)
      
      # make matrices, note that only predictor variables get normalized
      trainsetpvs <- as.matrix(trainsetpvs_normalized)
      testsetpvs <- as.matrix(testsetpvs_normalized)
      trainsetrv <- as.matrix(trainset$FLOW)
      testsetrv <- as.matrix(testset$FLOW)
      
      # make an array out of matrices
      trainset_array <- makearray(trainset, trainsetpvs)
      testset_array <- makearray(testset, testsetpvs)
      
      # reshape arrays into moving window for LSTM training
      trainset_array_window <- reshapearray(trainset_array, ncells)
      testset_array_window <- reshapearray(testset_array, ncells)
      
      # split into predictor variables and response variables
      trainsetpv_array <- trainset_array_window[ , , 1:nfeat]
      trainsetrv_array <- trainset_array_window[ , , nfeat+1]
      testsetpv_array <- testset_array_window[ , , 1:nfeat]
      testsetrv_array <- testset_array_window[ , , nfeat+1]
      
      # reset lr in case a callback in a previous loop has changed it
      lr <- 1e-3
      
      nnmodel %>%
        fit(trainsetpv_array, trainsetrv_array, epochs = nepochs, verbose=verbose, batch_size=320, validation_split = 0.2, shuffle = FALSE, callbacks = callbacks)
      
      predictions <- nnmodel %>% predict(testsetpv_array)
      predictions <- predictions[ , 1]
      results[[k]] <- cbind(obs=testsetrv_array[, ncells], pred=predictions)
      modgof[[k]] <- gof(as.data.frame(results[[k]])$pred, as.data.frame(results[[k]])$obs)
    }
    
    names_tbd <- rownames(modgof[[1]])
    modgof <- data.frame(matrix(unlist(modgof), nrow=length(modgof[[1]]), byrow=FALSE))
    rownames(modgof) <- names_tbd
  }
  
  if(groupingstyle=="logo"){
    list(results=results, gof=modgof)
  } else if(groupingstyle=="lmgo"){
    list(results=results, gof=modgof, kf=group)
  } else{
    print("not a valid groupingstyle!")
  }
}

