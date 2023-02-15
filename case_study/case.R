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
model_files <- c("forecast_ETAS_LM_FP32.dat.xz",
                 "forecast_ETES_FCM_FP32.dat.xz",
                 "forecast_STEP_LG_FP32.dat.xz",
                 "forecast-ensemble_SMA_FP32.dat.xz")
# Time stamps corresponding to model outputs
# (rows of the model output data)
time_stamps_file <- "meta_rows_dates.csv"
# Locations of grid cells corresponding to model outputs
# (columns of the model output data)
cell_file <- "meta_columns_cells.csv"
# Catalog of observed earthquakes
events_file <- "catalog.csv"

# Load necessary packages
library(Matrix)
library(maps)


# Load auxiliary functions
source(file.path(rpath, "case_functions.R"))

# Set model names and their colors
model_names <- c("LM", "FMC", "LG", "SMA")
model_colors <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
last_day <- list(DD = 20, MM = 5, YY = 2020)


###############
## Load data ##

# Load time stamps for the models
file_path <- file.path(dpath, time_stamps_file)
res <- load_times(file_path, last_day)
times <- res$times

# Load list of model forecasts
file_paths <- file.path(dpath, model_files)
models <- load_models(file_paths, res$tindex)
nmods <- length(models)

# Load the grid cell data (testing region)
file_path <- file.path(dpath, cell_file)
cells <- load_cells(file_path)

# Load data frame of M4+ events
file_path <- file.path(dpath, events_file)
events <- load_events(file_path, times)

# Filter the M4+ events for testing region
events <- filter_region(events, cells)

# Convert events data frame to observation matrix of the
# same format as the forecast matrices in models. Thus 
# any scoring function S can be applied to the components
# of these matrices. For example S( models[[i]], observations )
# gives a matrix of scores for all grid cells and days
ncells <- dim(cells)[1]
ndays <- dim(times)[1]
observations <- events2obs(events, ndays, ncells)


####################
## Plots and maps ##

# Compute daily and overall scores for the quadratic
# and the Poisson scoring function, see Section 5.2
scores_pois <- scores_quad <- matrix(0, nrow = ndays, ncol = nmods)
for (i in 1:nmods) {
  scores_pois[ ,i] <- rowSums( Spois(models[[i]], observations) )
  scores_quad[ ,i] <- rowSums( Squad(models[[i]], observations) )
}

# Create lateX table for the overall mean scores
file_path <- file.path(fpath, "score_table.tex")
scores2teX(scores_pois, scores_quad, model_names, file_path)

# Create plot of daily scores (Poisson)
file_path <- file.path(fpath, "plot_pois_time_log.pdf")
plot_scores(scores_pois, times, model_names, model_colors, file_path, events = events) 

# Create plot of daily scores (Quadratic)
file_path <- file.path(fpath, "plot_quad_time_log.pdf")
plot_scores(scores_quad, times, model_names, model_colors, file_path, events = events) 


# Plot maps of score differences
for (i in 2:nmods) {
  score_diff <- colMeans( Spois(models[[1]], observations)
                          - Spois(models[[i]], observations) )
  file_path <- file.path(fpath, paste0("map_score_diff_1", i, ".pdf"))
  plot_diffs_map(score_diff, cells, file_path)
}


# Plot ACF of score differences
file_path <- file.path(fpath, "plot_score_acf_pois.pdf")
plot_score_diffs_acf(models, observations, Spois, file_path)
file_path <- file.path(fpath, "plot_score_acf_quad.pdf")
plot_score_diffs_acf(models, observations, Squad, file_path)

# Create lateX tables with mean, variance, and detectable differences for all
# model combinations
pois_diffs <- calculate_means_and_vars(models, observations, Spois)
quad_diffs <- calculate_means_and_vars(models, observations, Squad)

file_path <- file.path(fpath, "table_sample_size_pois.tex")
sample_size_table(pois_diffs$means, pois_diffs$vars, pois_diffs$names, ndays,
                  "Poisson score", file_path)
file_path <- file.path(fpath, "table_sample_size_quad.tex")
sample_size_table(quad_diffs$means, quad_diffs$vars, quad_diffs$names,
                  ndays, "Quadratic score", file_path, scaling = 10)


# Plot map of aggregated score differences
# Compute neighborhood matrix
k <- 5
nmat <- neigh_mat(cells, k)

# Aggregate forecast models and observations
for (i in 1:nmods) {
  models[[i]] <- as.matrix( models[[i]] %*% nmat )
  gc()
}
observations <- observations %*% nmat

for (i in 2:nmods) {
  score_diff <- colMeans( Spois(models[[1]], observations)
                          - Spois(models[[i]], observations) )
  file_path <- file.path(fpath, paste0("map_score_diff_1", i, "_agg", k, ".pdf"))
  plot_diffs_map(score_diff, cells, file_path)
}

