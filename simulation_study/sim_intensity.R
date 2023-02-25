######################################################
## Simulation experiments for the intensity measure ##

#################
## Preparation ##

# Path for r scripts
rpath <- getwd()
  
# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/figures"
###
### Set to getwd() later

# Set a seed
set.seed(2020)

# Specify some global parameters
N <- 100   # sample size
M <- 500   # number of repetitions
alpha <- 0.05   # level for DM tests
table_color <- "cyan"

# Define forecast intensity functions
f0 <- function(x,y) 6 * sqrt(x^2 + y^2)
f1 <- function(x,y) 7.8 * sqrt((x-0.2)^2 + (y-0.1)^2)
f2 <- function(x,y) 2.3 * abs(x + 3 * y)
f3 <- function(x,y) 10 * sqrt((x-0.2)^2 + (y-0.1)^2)
f4 <- function(x,y) 7.5 * exp( - 3 * ( (x-0.6)^2 + (y-0.6)^2) )
f5 <- function(x,y) 2 * ( 1/sqrt(1.2 - x) + 2*(1-y) )
forecast_functions <- list(f0, f1, f2, f3, f4, f5)
n_forecast <- length(forecast_functions)

# Specify point processes
DPPscale <- 0.06
DPPvar <- f0(1,1)   # has to be maximum of f0
LGCPcov <- "exp"
LGCPvar <- 1/4
LGCPscale <- 1/5
LGCPmu <- function(x,y) log(f0(x,y)) - LGCPvar/2
TPPscale <- 0.05
TPPmu <- 1.5
TPPkappa <- function(x,y)  f0(x,y)/TPPmu
# Specify names of our four point process models
model_names <- c("Poisson", "DPP", "LGCP", "Thomas")

# Load necessary packages
library(spatstat)
library(cubature)
library(RandomFields)

# Load auxiliary functions
source(file.path(rpath, "sim_functions.R"))

# Compute integrals of forecast intensities
forecast_norms <- sapply(forecast_functions, integrate2D)
forecast_index <- 0:5


###########################
## Simulations and plots ##

# Do 2D plot of intensity forecasts
file_path <- file.path(fpath, "plot_intensities.pdf")
plotIntensities(forecast_functions, file_path)


## 1) Scoring function S1

# Save mean score differences of each experiment for boxplot
score_diffs1 <- list()

## 1.1) Poisson point process
res <- experiment("Poisson", ScoreIntensity1, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s1_Poisson.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
# Compute score differences for boxplots
score_diffs1[[1]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## 1.2) Determinantal point process
res <- experiment("DPP", ScoreIntensity1, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s1_DPP.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
# Compute score differences for boxplots
score_diffs1[[2]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## 1.3) Log-Gaussian Cox process
res <- experiment("LGCP", ScoreIntensity1, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s1_LGCP.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
# Compute score differences for boxplots
score_diffs1[[3]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## 1.4) Thomas cluster process
res <- experiment("Thomas", ScoreIntensity1, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s1_Thomas.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
# Compute score differences for boxplots
score_diffs1[[4]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## Create boxplot for scoring function S1
file_path <- file.path(fpath, "boxplot_intensity_s1.pdf")
diffmin <- sapply(score_diffs1, min)
diffmax <- sapply(score_diffs1, max)
ylim <- rep(list(c(min(diffmin), max(diffmax))), 4)
myBoxplot(score_diffs1, 1:5, model_names, ylim, file_path)


## 2) Scoring function S2

# Save mean score differences of each experiment for boxplot
score_diffs2 <- list()
# Save DM results of each experiment for the convergence experiments
# in part 3)
DM_results <- list()

## 2.1) Poisson point process
res <- experiment("Poisson", ScoreIntensity2, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s2_Poisson.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
DM_results[[1]] <- res$DM
# Compute score differences for boxplots
score_diffs2[[1]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## 2.2) Determinantal point process
res <- experiment("DPP", ScoreIntensity2, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s2_DPP.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
DM_results[[2]] <- res$DM
# Compute score differences for boxplots
score_diffs2[[2]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## 2.3) Log-Gaussian Cox process
res <- experiment("LGCP", ScoreIntensity2, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s2_LGCP.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
DM_results[[3]] <- res$DM
# Compute score differences for boxplots
score_diffs2[[3]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## 2.4) Thomas cluster process
res <- experiment("Thomas", ScoreIntensity2, forecast_functions, forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_s2_Thomas.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
DM_results[[4]] <- res$DM
# Compute score differences for boxplots
score_diffs2[[4]] <- res$Scores[ ,2:n_forecast] - res$Scores[ ,rep(1, n_forecast-1)]


## Create boxplot for scoring function S2
file_path <- file.path(fpath, "boxplot_intensity_s2.pdf")
diffmin <- sapply(score_diffs2, min)
diffmax <- sapply(score_diffs2, max)
ylim <- rep(list(c(min(diffmin), max(diffmax))), 4)
myBoxplot(score_diffs2, 1:5, model_names, ylim, file_path)


## 3) Grid cell approximation

# Define the number of sets in the axis partition
n_partition <- 2^(1:6)

# Compute grid cell forecasts for all values of n
forecast_cell_integrals <- list()
for (i in 1:n_forecast) {
  forecast_cell_integrals[[i]] <- lapply(n_partition, function(k) cellIntegrals(forecast_functions[[i]], k))
}

# Save the DM results for each experiment and all values of n
DM_cells <- list()
# Define the index of the forecast for which the convergence of the
# preferring probabilities should be studied. If =1 then calculate the
# probability of preferring f_0 over the other forecasts 
forecast_base <- 1
  
## 3.1) Poisson point process
DM_cells[[1]] <- experiment_cells("Poisson", forecast_cell_integrals, n_partition, forecast_base, M, N)

## 3.2) Determinantal point process
DM_cells[[2]] <- experiment_cells("DPP", forecast_cell_integrals, n_partition, forecast_base, M, N)

## 3.3) Log-Gaussian Cox process
DM_cells[[3]] <- experiment_cells("LGCP", forecast_cell_integrals, n_partition, forecast_base, M, N)

## 3.4) Thomas cluster process
DM_cells[[4]] <- experiment_cells("Thomas", forecast_cell_integrals, n_partition, forecast_base, M, N)


## Create plot for DM frequency convergence
file_path <- file.path(fpath, "plot_DM_convergence.pdf")
plotDMconv(DM_cells, DM_results, 1:6, forecast_base, model_names, file_path)
