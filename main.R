set.seed(3232019)

# libraries
library(reshape2)
library(keras)
library(dismo)

# data
moddf <- readRDS('Intermediary Data/moddf.rds')

# functions
source("functions.R")

# make an instance of keras and tensorflow
install_keras()

# # nn experiments 
# nns_resub_results <- nnscv(moddf, "resub")
# saveRDS(nns_resub_results, file="Output Data/rds/nns_resub_results.RDS")
# nns_2fold_results <- nnscv(moddf, 2L)
# saveRDS(nns_2fold_results, file="Output Data/rds/nns_2fold_results.RDS")
# nns_5fold_results <- nnscv(moddf, 5L)
# saveRDS(nns_5fold_results, file="Output Data/rds/nns_5fold_results.RDS")
# nns_10fold_results <- nnscv(moddf, 10L)
# saveRDS(nns_10fold_results, file="Output Data/rds/nns_10fold_results.RDS")
# nns_logo_results <- nnscv(moddf, "logo")
# saveRDS(nns_logo_results, file="Output Data/rds/nns_logo_results.RDS")
# nns_lmgo_results <- nnscv(moddf, "lmgo")
# saveRDS(nns_lmgo_results, file="Output Data/rds/nns_lmgo_results.RDS")

# # save nn experiments
# nns_resub_results <- readRDS("Output Data/rds/nns_resub_results.RDS")
# nns_2fold_results <- readRDS("Output Data/rds/nns_2fold_results.RDS")
# nns_5fold_results <- readRDS("Output Data/rds/nns_5fold_results.RDS")
# nns_10fold_results <- readRDS("Output Data/rds/nns_10fold_results.RDS")
# nns_logo_results <- readRDS("Output Data/rds/nns_logo_results.RDS")
# nns_lmgo_results <- readRDS("Output Data/rds/nns_lmgo_results.RDS")

# # lstm experiments 
lstm_logo_exp1_results <- lstmcv(moddf, "logo")
# lstm_logo_exp2_results <- lstmcv(moddf, "logo", lstm_units = 256)
# lstm_logo_exp3_results <- lstmcv(moddf, "logo", lstm_units = 256, ncells=50)
# lstm_logo_exp4_results <- lstmcv(moddf, "logo", lstm_units = 256, ncells=50, dropout_rate = 0.2)
# lstm_logo_exp5_results <- lstmcv(moddf, "logo", lstm_units = 256, ncells=50, dropout_rate = 0, factor=1e-4)
# lstm_logo_exp6_results <- lstmcv(moddf, "logo", lstm_units = 256, ncells=50, dropout_rate = 0, factor=NULL, patience_early = 100)
# lstm_logo_exp7_results <- lstmcv(moddf, "logo", lstm_units = 256, ncells=50, dropout_rate = 0, factor=NULL, patience_early = NULL)
# 
# lstm_lmgo_exp1_results <- lstmcv(moddf, "lmgo")
# lstm_lmgo_exp2_results <- lstmcv(moddf, "lmgo", lstm_units = 256)
# lstm_lmgo_exp3_results <- lstmcv(moddf, "lmgo", lstm_units = 256, ncells=50)
# lstm_lmgo_exp4_results <- lstmcv(moddf, "lmgo", lstm_units = 256, ncells=50, dropout_rate = 0.2)
# lstm_lmgo_exp5_results <- lstmcv(moddf, "lmgo", lstm_units = 256, ncells=50, dropout_rate = 0, factor=1e-4)
# lstm_lmgo_exp6_results <- lstmcv(moddf, "lmgo", lstm_units = 256, ncells=50, dropout_rate = 0, factor=NULL, patience_early = 100)
# lstm_lmgo_exp7_results <- lstmcv(moddf, "lmgo", lstm_units = 256, ncells=50, dropout_rate = 0, factor=NULL, patience_early = NULL)
# 
# # save experiments
saveRDS(lstm_logo_exp1_results, file="Output Data/rds/lstm_logo_exp1_results.RDS")
# saveRDS(lstm_logo_exp2_results, file="Output Data/rds/lstm_logo_exp2_results.RDS")
# saveRDS(lstm_logo_exp3_results, file="Output Data/rds/lstm_logo_exp3_results.RDS")
# saveRDS(lstm_logo_exp4_results, file="Output Data/rds/lstm_logo_exp4_results.RDS")
# saveRDS(lstm_logo_exp5_results, file="Output Data/rds/lstm_logo_exp5_results.RDS")
# saveRDS(lstm_logo_exp6_results, file="Output Data/rds/lstm_logo_exp6_results.RDS")
# saveRDS(lstm_logo_exp7_results, file="Output Data/rds/lstm_logo_exp7_results.RDS")
# saveRDS(lstm_logo_exp8_results, file="Output Data/rds/lstm_logo_exp8_results.RDS")
# 
# saveRDS(lstm_lmgo_exp1_results, file="Output Data/rds/lstm_lmgo_exp1_results.RDS")
# saveRDS(lstm_lmgo_exp2_results, file="Output Data/rds/lstm_lmgo_exp2_results.RDS")
# saveRDS(lstm_lmgo_exp3_results, file="Output Data/rds/lstm_lmgo_exp3_results.RDS")
# saveRDS(lstm_lmgo_exp4_results, file="Output Data/rds/lstm_lmgo_exp4_results.RDS")
# saveRDS(lstm_lmgo_exp5_results, file="Output Data/rds/lstm_lmgo_exp5_results.RDS")
# saveRDS(lstm_lmgo_exp6_results, file="Output Data/rds/lstm_lmgo_exp6_results.RDS")
# saveRDS(lstm_lmgo_exp7_results, file="Output Data/rds/lstm_lmgo_exp7_results.RDS")
# saveRDS(lstm_lmgo_exp8_results, file="Output Data/rds/lstm_lmgo_exp8_results.RDS")
