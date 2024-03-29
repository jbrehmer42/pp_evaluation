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
table_color <- "cyan"

# Specify point processes
LGCPcov <- "gauss"
LGCPvar <- log(2)
LGCPscale <- 1/20
LGCPmu <- log(lambda) - LGCPvar/2
DPPscale <- 0.06
# Specify names of our three point process models
model_names <- c("LGCP", "Poisson", "DPP")

# Define forecast product densities via functions on [0,\infty)
# These functions take the role of \rho_0 in Example B.6
f1 <- function(r) exp( 2 * LGCPmu + LGCPvar * (1 + exp(- (r/LGCPscale)^2)) )
f2 <- function(r) exp( 2 * LGCPmu + LGCPvar * (1 + exp(- r/LGCPscale)) )
f3 <- function(r) rep(lambda^2, length(r))
f4 <- function(r) lambda^2 * (1 - exp(-2 * r/DPPscale) )
f5 <- function(r) lambda^2 * (1 - exp(-2 * (r/DPPscale)^2) )
forecast_functions <- list(f1, f2, f3, f4, f5)
n_forecast <- length(forecast_functions)

# Load necessary packages
library(spatstat)
library(cubature)
library(RandomFields)

# Load auxiliary functions
source(file.path(rpath, "sim_functions.R"))

# Compute integrals of forecast of the product density
# which corresponds to the functions f1, ..., f5
forecast_norms <- sapply(forecast_functions, integrate4D)
forecast_index <- 1:5


###########################
## Simulations and plots ##

# Plot of product densities
file_path <- file.path(fpath, "plot_product_densities.pdf")
plot_product_densities(forecast_functions, file_path)

# Save mean scores of each experiment for boxplot
scores <- list()


## 1) Log-Gaussian Cox process
res <- experiment("homLGCP", S_product_density, forecast_functions,
                  forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_prod_LGCP.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
scores[[1]] <- res$Scores


## 2) Poisson point process
res <- experiment("homPoisson", S_product_density, forecast_functions,
                  forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_prod_Poisson.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
scores[[2]] <- res$Scores


## 3) Determinantal point process
res <- experiment("homDPP", S_product_density, forecast_functions,
                  forecast_norms, M, N)
# Print the DM table
file_path <- file.path(fpath, "DM_table_prod_DPP.tex")
array2teX(res$DM, forecast_index, file_path, color = table_color)
scores[[3]] <- res$Scores


## Create boxplot
file_path <- file.path(fpath, "boxplot_product.pdf")
ylim <- c( min(sapply(scores, min)), max(sapply(scores, max)) )
ylim <- rep(list(ylim), 3)
forecast_names <- rep(list(paste0("f", 1:5)), 3)
for (i in 1:3) {
  forecast_names[[i]][2*(i-1)+1] <- paste0("f", 2*(i-1)+1, " (true)")
}
boxplot_score(scores, forecast_names , model_names, ylim, file_path,
              zero_line = FALSE)
