# Predicting Ungauged Basins for Daily Unimpaired Flows
A record of data gathering, processing, modelling with differnet loss functions, and resampling methods for the predicting ungauged basins (PUB) problem. All code was written in R and uses many libraries (see citations in \*.Rmd files).

## Installation
Install R and Rstudio. Clone or donwload this repository. Open Rstudio and create a R Project in the cloned directory. Open the .Rmd files in the R Project environment. Install all missing libraries with install.packages() in the Rstudio console, e.g.:  
```bash
install.packages("ratser")
```

No GPU is required for running keras and training neural networks.

## Usage 
All code used to process data and for modelling is in \*.Rmd files. 

You will need an "Input Data" folder to run "ch1 data prep.Rmd." This folder was too large to host on GitHub. Email me at white.elaheh@gmail.com for these files. 

The "Intermediary Data" folder was too large to host on Github, and includes the processed data coming from "ch1 data prep.Rmd." This is the dataframe used in all the models: "moddf.rds."

The "Libraries" folder includes r code used to calculate model measures of fit. I modified the bR2 function of the HydroGOF library, which is why I am loading the functions rather than loading the library. 

The "Output Data" folder contains all plots from the code but is missing a "rds" folder too large to host on GitHub. "Output Data/rds" saved all the models produced by the code in "ch2 visualizations.Rmd", "ch3 loss functions.Rmd", and "ch4 resampling.Rmd."

Once moddf.rds is made from running "ch1 data prep.Rmd", other \*.Rmd files can be ran in what ever order is wanted. 

"ml overview.pptx" is a slideshow explaining the concepts behind modeling for PUB. 

## Support 
Email me at white.elaheh@gmail.com. 

## Roadmap 
Future work can improve predictions by: 
1) including more basins in the study to give more diversity in basin characteristics.
2) adding more predictor variables to the data: e.g., min and max temperature, solar radiation, and vapor pressure, average seasonal precipitations, average annual temperatures, precent coniferous or coniferous/deciduous land cover, cumulative monthly (or weekly) precipitation.  
3) changing the neural network architecture to include many basin observations as inputs (although this may be tricky with predicting one basin at a time in the LOGO cross validation method) and include a LSTM to take into account the sequential nature of the observations. 
4) including network information with a more sophisticated model suited to handling network flow data (e.g., traffic prediciton models). 

## Project Status 
Development has stopped as of 11/30/2020. 