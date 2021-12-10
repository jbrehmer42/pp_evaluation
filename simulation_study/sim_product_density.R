####################################################
## Simulation experiments for the product density ##

#################
## Preparation ##

# Path for r scripts
rpath <- getwd()

# Path for figures
fpath <- getwd()

# Set a seed
set.seed(2020)

# Specify some global parameters
N <- 30    # sample size
M <- 500   # number of repetitions
alpha <- 0.05   # level for DM tests
lambda <- 25
tabColor <- "cyan"

# Specify point processes
LGCPcov <- "gauss"
LGCPvar <- log(2)
LGCPscale <- 1/20
LGCPmu <- log(lambda) - LGCPvar/2
DPPscale <- 0.06
# Specify names of our three point process models
mnames <- c("LGCP", "Poisson", "DPP")

# Define forecast product densities via functions on [0,\infty)
# These functions take the role of \rho_0 in Example B.6
f1 <- function(r) exp( 2 * LGCPmu + LGCPvar * (1 + exp(- (r/LGCPscale)^2)) )
f2 <- function(r) exp( 2 * LGCPmu + LGCPvar * (1 + exp(- r/LGCPscale)) )
f3 <- function(r) rep(lambda^2, length(r))
f4 <- function(r) lambda^2 * (1 - exp(-2 * r/DPPscale) )
f5 <- function(r) lambda^2 * (1 - exp(-2 * (r/DPPscale)^2) )
fcast_funs <- list(f1, f2, f3, f4, f5)
nfcast <- length(fcast_funs)

# Load necessary packages
library(spatstat)
library(cubature)
library(RandomFields)

# Load auxiliary functions
source(paste0(rpath, "/sim_functions.R"))

# Compute integrals of forecast of the product density
# which corresponds to the functions f1, ..., f5
fcast_norm <- sapply(fcast_funs, integrate4D)
fcast <- 1:5


###########################
## Simulations and plots ##

# Plot of product densities
filePath <- paste(fpath, "plot_product_densities.pdf", sep = "/")
plotProduct(fcast_funs, filePath)

# Save mean scores of each experiment for boxplot
scorelist <- list()


## 1) Log-Gaussian Cox process
res <- experiment("homLGCP", ScoreProductDensity, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- paste(fpath, "DM_table_prod_LGCP.tex", sep = "/")
array2teX(res$DM, fcast, filePath, color = tabColor)
scorelist[[1]] <- res$Scores


## 2) Poisson point process
res <- experiment("homPoisson", ScoreProductDensity, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- paste(fpath, "DM_table_prod_Poisson.tex", sep = "/")
array2teX(res$DM, fcast, filePath, color = tabColor)
scorelist[[2]] <- res$Scores


## 3) Determinantal point process
res <- experiment("homDPP", ScoreProductDensity, fcast_funs, fcast_norm, M, N)
# Print the DM table
filePath <- paste(fpath, "DM_table_prod_DPP.tex", sep = "/")
array2teX(res$DM, fcast, filePath, color = tabColor)
scorelist[[3]] <- res$Scores


## Create boxplot
filePath <- paste(fpath, "boxplot_product.pdf", sep = "/")
ylim <- c( min(sapply(scorelist, min)), max(sapply(scorelist, max)) )
ylim <- rep(list(ylim), 3)
fnames <- rep(list(paste0("f", 1:5)), 3)
for (i in 1:3) fnames[[i]][2*(i-1)+1] <- paste0("f", 2*(i-1)+1, " (true)")
myBoxplot(scorelist, fnames , mnames, ylim, filePath, zline = F)

