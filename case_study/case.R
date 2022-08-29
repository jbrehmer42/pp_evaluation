#################
## Preparation ##

# Path for r scripts
rpath <- getwd()

# Path for figures
fpath <- "/media/myData/Plots"

# Path for data
# The files containing the forecasts and observations should
# be located in this folder
dpath <- "/media/myData/EQData"
###
### Set to getwd() later

# Set file names (default names)
# The forecast model outputs are arrays with time in rows
# and grid cells in the columns
modelNames <- c("ETAS_LM.txt.xz", "ETES_FMC.txt.xz",
                "STEP_LG.txt.xz", "Bayesian_corr_27_10.txt.xz")
# Time stamps corresponding to model outputs
# (rows of the model output data)
timestampName <- "meta_rows.txt"
# Locations of grid cells corresponding to model outputs
# (columns of the model output data)
cellName <- "meta_column.csv"
# Catalog of observed earthquakes
eventsName <- "meta_catalogo.txt"

# Load necessary packages
library(Matrix)
library(maps)


# Load auxiliary functions
source(file.path(rpath, "case_functions.R"))

# Set model names and their colors
mnames <- c("LM", "FMC", "LG", "SMA")
mcols <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
lastday <- list(DD = 20, MM = 5, YY = 2020)


###############
## Load data ##

# Load time stamps for the models
filePath <- file.path(dpath, timestampName)
res <- load_times(filePath, lastday)
times <- res$times

# Load list of model forecasts
filePaths <- file.path(dpath, modelNames)
models <- load_models(filePaths, res$tindex)
nmods <- length(models)

# Load the grid cell data (testing region)
filePath <- file.path(dpath, cellName)
cells <- load_cells(filePath)

# Load data frame of M4+ events
filePath <- file.path(dpath, eventsName)
events <- load_events(filePath, times)

# Filter the M4+ events for testing region
events <- filterRegion(events, cells)

# Convert events data frame to observation matrix of the
# same format as the forecast matrices in models. Thus 
# any scoring function S can be applied to the components
# of these matrices. For example S( models[[i]], obs )
# gives a matrix of scores for all grid cells and days
ncells <- dim(cells)[1]
ndays <- dim(times)[1]
obs <- events2obs(events, ndays, ncells)


####################
## Plots and maps ##

# Compute daily and overall scores for the quadratic
# and the Poisson scoring function, see Section 5.2
scores_pois <- scores_quad <- matrix(0, nrow = ndays, ncol = nmods)
for (i in 1:nmods) {
  scores_pois[ ,i] <- rowSums( Spois(models[[i]], obs) )
  scores_quad[ ,i] <- rowSums( Squad(models[[i]], obs) )
}

# Create lateX table for the overall mean scores
filePath <- file.path(fpath, "score_table.tex")
scores2teX(scores_pois, scores_quad, mnames, filePath)

# Create plot of daily scores (Poisson)
filePath <- file.path(fpath, "plot_pois_time_log.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events = events) 

# Create plot of daily scores (Quadratic)
filePath <- file.path(fpath, "plot_quad_time_log.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events = events) 


# Plot maps of score differences
for (i in 2:nmods) {
  sdiff <- colMeans( Spois(models[[1]], obs) - Spois(models[[i]], obs) )
  filePath <- file.path(fpath, paste0("map_score_diff_1", i, ".pdf"))
  diffMap(sdiff, cells, filePath)
}

## Plot map of aggregated score differences
# Compute neighborhood matrix
k <- 5
nmat <- neigh_mat(cells, k)

# Aggregate forecast models and observations
for (i in 1:nmods) {
  models[[i]] <- as.matrix( models[[i]] %*% nmat )
  gc()
}
obs <- obs %*% nmat

for (i in 2:nmods) {
  sdiff <- colMeans( Spois(models[[1]], obs) - Spois(models[[i]], obs) )
  filePath <- file.path(fpath, paste0("map_score_diff_1", i, "_agg", k, ".pdf"))
  diffMap(sdiff, cells, filePath)
}

