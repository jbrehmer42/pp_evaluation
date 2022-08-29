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
tabColor <- "cyan"

# Define forecast intensity functions
f0 <- function(x,y) 6 * sqrt(x^2 + y^2)
f1 <- function(x,y) 7.8 * sqrt((x-0.2)^2 + (y-0.1)^2)
f2 <- function(x,y) 2.3 * abs(x + 3 * y)
f3 <- function(x,y) 10 * sqrt((x-0.2)^2 + (y-0.1)^2)
f4 <- function(x,y) 7.5 * exp( - 3 * ( (x-0.6)^2 + (y-0.6)^2) )
f5 <- function(x,y) 2 * ( 1/sqrt(1.2 - x) + 2*(1-y) )
fcast_funs <- list(f0, f1, f2, f3, f4, f5)
nfcast <- length(fcast_funs)

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
mnames <- c("Poisson", "DPP", "LGCP", "Thomas")

# Load necessary packages
library(spatstat)
library(cubature)
library(RandomFields)

# Load auxiliary functions
source(file.path(rpath, "sim_functions.R"))

# Compute integrals of forecast intensities
fcast_norm <- sapply(fcast_funs, integrate2D)
fcast <- 0:5


###########################
## Simulations and plots ##

# Do 2D plot of intensity forecasts
filePath <- file.path(fpath, "plot_intensities.pdf")
plotIntensities(fcast_funs, filePath)


## 1) Scoring function S1

# Save mean score differences of each experiment for boxplot
difflist1 <- list()

## 1.1) Poisson point process
res <- experiment("Poisson", ScoreIntensity1, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s1_Poisson.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
# Compute score differences for boxplots
difflist1[[1]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## 1.2) Determinantal point process
res <- experiment("DPP", ScoreIntensity1, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s1_DPP.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
# Compute score differences for boxplots
difflist1[[2]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## 1.3) Log-Gaussian Cox process
res <- experiment("LGCP", ScoreIntensity1, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s1_LGCP.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
# Compute score differences for boxplots
difflist1[[3]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## 1.4) Thomas cluster process
res <- experiment("Thomas", ScoreIntensity1, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s1_Thomas.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
# Compute score differences for boxplots
difflist1[[4]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## Create boxplot for scoring function S1
filePath <- file.path(fpath, "boxplot_intensity_s1.pdf")
diffmin <- sapply(difflist1, min)
diffmax <- sapply(difflist1, max)
ylim <- rep(list(c(min(diffmin), max(diffmax))), 4)
myBoxplot(difflist1, 1:5, mnames, ylim, filePath)


## 2) Scoring function S2

# Save mean score differences of each experiment for boxplot
difflist2 <- list()
# Save DM results of each experiment for the convergence experiments
# in part 3)
DMlist <- list()

## 2.1) Poisson point process
res <- experiment("Poisson", ScoreIntensity2, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s2_Poisson.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
DMlist[[1]] <- res$DM
# Compute score differences for boxplots
difflist2[[1]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## 2.2) Determinantal point process
res <- experiment("DPP", ScoreIntensity2, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s2_DPP.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
DMlist[[2]] <- res$DM
# Compute score differences for boxplots
difflist2[[2]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## 2.3) Log-Gaussian Cox process
res <- experiment("LGCP", ScoreIntensity2, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s2_LGCP.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
DMlist[[3]] <- res$DM
# Compute score differences for boxplots
difflist2[[3]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## 2.4) Thomas cluster process
res <- experiment("Thomas", ScoreIntensity2, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- file.path(fpath, "DM_table_s2_Thomas.tex")
array2teX(res$DM, fcast, filePath, color = tabColor)
DMlist[[4]] <- res$DM
# Compute score differences for boxplots
difflist2[[4]] <- res$Scores[ ,2:nfcast] - res$Scores[ ,rep(1, nfcast-1)]


## Create boxplot for scoring function S2
filePath <- file.path(fpath, "boxplot_intensity_s2.pdf")
diffmin <- sapply(difflist2, min)
diffmax <- sapply(difflist2, max)
ylim <- rep(list(c(min(diffmin), max(diffmax))), 4)
myBoxplot(difflist2, 1:5, mnames, ylim, filePath)


## 3) Grid cell approximation

# Define the number of sets in the axis partition
nset <- 2^(1:6)

# Compute grid cell forecasts for all values of n
fcast_cell <- list()
for (i in 1:nfcast) fcast_cell[[i]] <- lapply(nset, function(k) cellIntegrals(fcast_funs[[i]], k))

# Save the DM results for each experiment and all values of n
approxlist <- list()
# Define the index of the forecast for which the convergence of the
# preferring probabilities should be studied. If =1 then calculate the
# probability of preferring f_0 over the other forecasts 
fbase <- 1
  
## 3.1) Poisson point process
approxlist[[1]] <- experiment_cells("Poisson", fcast_cell, nset, fbase, M, N)

## 3.2) Determinantal point process
approxlist[[2]] <- experiment_cells("DPP", fcast_cell, nset, fbase, M, N)

## 3.3) Log-Gaussian Cox process
approxlist[[3]] <- experiment_cells("LGCP", fcast_cell, nset, fbase, M, N)

## 3.4) Thomas cluster process
approxlist[[4]] <- experiment_cells("Thomas", fcast_cell, nset, fbase, M, N)


## Create plot for DM frequency convergence
filePath <- file.path(fpath, "plot_DM_convergence.pdf")
plotDMconv(approxlist, DMlist, 1:6, fbase, mnames, filePath)
