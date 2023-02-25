######################################
## Case study - Auxiliary functions ##

# Poisson scoring function
# Caution: Have to ensure that reports x are never zero
Spois <- function(x, y)  -y * log(x) + x

# Quadratic scoring function
Squad <- function(x,y) (y - x)^2

# Detectable difference estimate of Lehr (1992)
lehr_estimate <- function(var, n) sqrt(8 * var / n)


# Load time stamps of the model outputs
load_times <- function(file_path, last_day) {
  # Input values:
  # file_path - File path for the time stamps file
  # last_day  - Last day for which forecasts can be
  #             evaluated.
  # Read time stamps from file
  raw_strings <- read.table(file_path, skip = 1, sep = ";")
  # Convert time stamps to data frame of integers
  times <- sapply(raw_strings,
                  function (x) unlist(regmatches(x, gregexpr("(\\d+)", x))))
  times <- as.data.frame(matrix(as.integer(times), ncol = 6, byrow = TRUE))
  names(times) <-  c("YY", "MM", "DD", "H", "M", "S")
  # Filter out days with multiple model runs: Compute
  # index of days with only one model run. Only the
  # first time stamp of the model runs is retained
  n <- dim(times)[1]
  tindex <- rep(T, n)
  for (i in 2:n) {
    # Check whether previous day agrees with current day
    tindex[i] <- !all(times[i-1, 1:3] == times[i, 1:3])
  }
  # Adjust for last day for which data is available. All
  # days after this day will be deleted
  match_last_day <- (times$DD == last_day$DD) & (times$MM == last_day$MM) &
    (times$YY == last_day$YY)
  if (any(match_last_day)) {
    tindex[max(which(match_last_day)+1):length(tindex)] <- FALSE
  }
  times <- times[tindex, ]
  # Return new time stamps and index of corresponding rows
  return(list(times = times, tindex = tindex))
}


# Load forecast values of the models
load_models <- function(file_paths, tindex) {
  # Input values:
  # file_paths - File paths for the model outputs
  # tindex     - Index of model runs (rows) to be used
  # Assume that the file tindex applies to all
  # forecast model files
  # Read forecast values from files
  nmod <- length(file_paths)
  models <- list()
  for (i in 1:nmod) {
    # Rows are days
    # Columns are grid cells
    models[[i]] <- as.matrix( read.table( xzfile(file_paths[i]) ) )
    attr(models[[i]], "dimnames") <- NULL
    # delete rows corresponding to multiple model runs on one day
    models[[i]] <- models[[i]][tindex, ]
  }
  # Return list of model output matrices
  return(models)
}


# Load grid cells (testing region)
load_cells <- function(file_path) {
  # Input values:
  # file_path - File path for the grid cells file
  # Read grid cell file
  cells <- read.csv(file_path, header = FALSE, skip = 1,
                    col.names = c("LON", "LAT"))
  # Start numbering the grid cells with 1 (just for
  # convenience and interpretability)
  cells$N <- 1:dim(cells)[1]
  # Add x-y-coordinates for cells
  xvals <- sort(unique(cells$LON))
  yvals <- sort(unique(cells$LAT))
  cells$X <- sapply(cells$LON, function(x) which(x == xvals))
  cells$Y <- sapply(cells$LAT, function(x) which(x == yvals))
  return(cells)
}


# Load events data.frame
load_events <- function(file_path, times) {
  # Input values:
  # file_path - File path for the events file
  # times     - Time stamps of testing period
  # Read events file
  events <- read.csv(file_path, header = FALSE, skip = 1,
                     col.names = c("time", "LAT", "LON", "DEP", "MAG"))
  # Convert time stamps to data frame of integers
  time_column <- sapply(events[1],
                      function(x) unlist(regmatches(x, gregexpr("(\\d+)", x))))
  time_df <- as.data.frame(matrix(as.numeric(time_column), ncol = 7, byrow = T))
  names(time_df) <- c("YY", "MM", "DD", "H", "M", "S", "MS")
  events <- cbind(time_df[1:6], events[-1])
  # Use only events with magnitude M >= 4
  # Using M >= 4 is equivalent to using M >= 3.95
  M4ind <- (events$MAG >= 4)
  events <- events[M4ind, ]
  # Assign the day number (time index TI) to events. Days
  # of the testing period (times) are consecutively numbered
  n <- dim(events)[1]
  time_index <- rep(-1, n)
  for (i in 1:n) {
    DD <- events$DD[i]
    MM <- events$MM[i]
    YY <- events$YY[i]
    gind <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY)
    if (any(gind)) time_index[i] <- which(gind)
  }
  # Erase events which do not occur during the testing period
  events <- events[(time_index > 0), ]
  # Add column time index (TI) to the events data frame
  events$TI <- time_index[time_index > 0]
  return(events)
}


# Filter out events which do not fall into the 
# testing region. Assign a grid cell number to
# the remaining events
filter_region <- function(events, cells) {
  # Input values:
  # events - Data frame of observed events
  # cells  - Data frame of grid cells
  # Specify edge length of grid cells
  size_LON <- 0.1
  size_LAT <- 0.1
  # Assign cell number to every event. Cell numbers
  # start with 1. If an event does not fall into the
  # testing region its cell number is -1 (and it is
  # filtered out)
  n <- dim(events)[1]
  ind <- rep(0, n)
  cell_ri <- cells$LON + 0.5 * size_LON
  cell_le <- cells$LON - 0.5 * size_LON
  cell_lo <- cells$LAT - 0.5 * size_LAT
  cell_up <- cells$LAT + 0.5 * size_LAT
  for (i in 1:n) {
    is_LON <- (cell_le < events$LON[i]) & (events$LON[i] <= cell_ri)
    is_LAT <- (cell_lo < events$LAT[i]) & (events$LAT[i] <= cell_up)
    if (any(is_LON & is_LAT) ) {
      # Assign cell number inside testing region
      ind[i] <- cells$N[is_LON & is_LAT]
    } else {
      # Set to -1 outside testing region
      ind[i] <- -1
    }
  }
  # Filter events and add cell number N
  events <- events[(ind > 0), ]
  events$N <- ind[ind > 0]
  return(events)
}


# Create an observation matrix from the events
# data frame
events2obs <- function(events, ndays, ncells) {
  # Input values:
  # events - Data frame of observed events
  # ndays  - Number of days in testing period
  # ncells - Number of cells in testing region
  # Create an observation matrix which can be directly
  # compared to the forecast model output matrices
  # Rows are days
  # Columns are grid cells
  observations <- Matrix(0, ncol = ncells, nrow = ndays, sparse = T)
  for (i in 1:ndays) {
    # Collect events in a 7-day period
    ind <- (events$TI >= i) & (events$TI < i + 7)
    if (any(ind)) {
      observations[i, ] <- tabulate(events$N[ind], nbins = ncells)
    } else next
  }
  return(observations)
}


# Print table of mean scores
scores2teX <- function(scores_pois, scores_quad, model_names, file_path) {
  # Input values:
  # scores_pois   - Matrix of daily scores
  # scores_quad   - Matrix of daily scores
  # model_names   - Names of the forecast models
  # file_path     - File path for the .tex file
  # Take score values and print teX code to
  # file to get values in tabular environment
  # Compute mean scores
  scores_pois <- colMeans(scores_pois)
  scores_quad <- colMeans(scores_quad)
  k <- length(scores_pois)
  begin <- "\\begin{tabular}{l cc}"
  head <- paste("Model", "pois", "quad \\\\", sep = " & ")
  # Write teX code to file
  write(begin, file_path)
  write("\\hline \\hline", file_path, append = T)
  write(head, file_path, append = T)
  write("\\hline", file_path, append = T)
  for (i in 1:k) {
    # Write score values
    row_pois <- sprintf('%.2f', scores_pois[i])
    row_quad <- sprintf('%.4f', scores_quad[i])
    row <- paste(model_names[i], row_pois, row_quad, sep = " & ")
    row <- paste0(row, " \\\\")
    write(row, file_path, append = T)
  }
  write("\\hline", file_path, append = T)
  write("\\end{tabular}", file_path, append = T)
}


# Plot daily scores of the models
plot_scores <- function(scores, times, model_names, model_colors, file_path,
                        events = NULL, logscale = T) {
  # Input values:
  # scores        - Matrix of daily scores
  # times         - Time stamps of the scores
  # model_names   - Model names for the legend
  # model_colors  - Array of colors for the lines
  # file_path     - File path for .pdf file
  # events        - Data frame of events (optional)
  # logscale      - Should scores be on log scale?
  ndays <- dim(times)[1]
  ylab <- "score"
  if(logscale) {
    scores <- log(scores)
    ylab <- paste(ylab, "(log scale)")
  }
  # Determine ylim
  ylim <- range(scores)
  # Ajust ylim if events are added to the plot
  if (hasArg(events)) {
    # Determine event circle coordinates
    evt_days <- unique(events$TI)
    evt_y <- ylim[1] - 3/96 * diff(ylim)
    ylim[1] <- ylim[1] - 4/96 * diff(ylim)
  }
  # Determine x-axis ticks and labels
  january1s <- (times$DD == 1) & (times$MM == 1)
  xlabels <- times$YY[january1s]
  xats <- which(january1s)
  pdf(file_path, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:ndays, rep(0, ndays), ylim = ylim, xlab = "days", col = "white",
       ylab = ylab, main = "", xaxt="n")
  for (i in 1:ncol(scores)) {
    lines(1:ndays, scores[ ,i], col = model_colors[i])
  }
  # Add event circles if possible
  if (hasArg(events))  points(evt_days, rep(evt_y, length(evt_days)), cex = 0.8)
  axis(1, at = xats, labels = xlabels)
  legend(-10, ylim[2], model_names, col = model_colors, lwd = 2)
  dev.off()
}


# Plot a map of score differences by individually
# coloring the grid cells
plot_diffs_map <- function(vals, cells, file_path, events = NULL, borders = T) {
  # Input values:
  # vals      - Numeric values corresponding to
  #            the grid cells
  # cells     - Data frame of grid cells
  # file_path - File path for .pdf file
  # events    - Data frame of events (optional)
  # borders   - Plot national borders? Requires 
  #             R-package "maps"
  # Set graphical paramters e.g. colors
  ylen <- 5
  nticks <- 7
  neg_col <- 0.66       # spec via hue in hsv colors
  pos_col <- 0          # spec via hue in hsv colors
  border_col <- rgb(0, 0, 0, alpha = 0.4) 
  # Transform values in vals into colors: Scale them to values in the
  # interval [-1, 1] and interpret these values as saturation
  val_abs <- 1.01 * max(abs(vals))
  scl <- (vals + val_abs) / (2 * val_abs)
  scl <- 2 * scl - 1
  col_scl <- rep("", length(vals))
  col_scl[scl >= 0] <- hsv(pos_col, scl[scl >= 0])
  col_scl[scl < 0]  <- hsv(neg_col, -scl[scl < 0])
  # Create pdf file with two parts
  pdf(file_path, width = 5, height = 4.4)
  layout(matrix(1:2, 1, 2, byrow = T), widths = c(0.78, 0.22))
  # Plot the map
  par(mar = 2/3 * rep(1,4))
  plot(1, 1, xlim = range(cells$LON), ylim = range(cells$LAT), col = "white",
       asp = 1.3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
  points(cells$LON, cells$LAT, pch = 15, col = col_scl, cex = 0.4)
  # Add national borders if wanted
  if (borders) map("world", fill = F, add = T, col = border_col)
  # Add locations of events if possible
  if (!is.null(events)) {
    binN <- unique(events$N)
    ind <- apply(as.matrix(cells$N), 1, function(x) any(x == binN))
    points(cells$LON[ind], cells$LAT[ind], pch = 5, col = border_col, 
           lwd = 0.5, cex = 0.4)
  }
  # Add color bar to the right of the map
  par(mar = c(1/3,0,1/3,2), mgp = c(-1,-2.4,-3.4), cex = 0.6)
  kk <- ylen * 20
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(-ylen, ylen), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labs <- sprintf("%.1e", seq(- val_abs, val_abs, len = nticks))
  labs[(nticks + 1)/2] <- "0"
  axis(4, at = seq(-ylen, ylen, len = nticks), labels = labs, las = 1)
  for (l in 1:(kk-1)) {
    sat <- l/(kk-1)
    rect(0, (l-1) * ylen/(kk - 1), 1/2,  l * ylen/(kk - 1), 
         col = hsv(pos_col, sat), border = NA)
    rect(0, - (l-1) * ylen/(kk - 1), 1/2,  -l * ylen/(kk - 1), 
         col = hsv(neg_col, sat), border = NA)
  }
  dev.off()
}


# Compute neighborhood matrix for spatial aggregation
# (See Section D in the Supplement)
neigh_mat <- function(cells, k) {
  # Input values:
  # cells - Data frame of grid cells
  # k     - Integer specifying the size of
  #         the neighborhood for aggregation
  ncells <- dim(cells)[1]
  # Aggregation will usually be done for small values
  # of k so "sparse = T" makes sense in most cases
  mat <- Matrix(0, nrow = ncells, ncol = ncells, sparse = T)
  for (i in 1:ncells) {
    x <- cells$X[i]
    y <- cells$Y[i]
    xind <- (cells$X <= x + k) & (cells$X >= x - k)
    yind <- (cells$Y <= y + k) & (cells$Y >= y - k)
    mat[i, ] <- as.numeric(xind & yind)
  }
  return(mat)
}


# Calculate mean and variance of the score differences for the computation
# of the sample size table in sample_size_table
calculate_means_and_vars <- function(models, observations, scoring_function) {
  # Input values:
  # models            - List of model outputs
  #                     matrices
  # observations      - Observation matrix
  # scoring_function  - Scoring function to compute
  #                     score differences
  means <- vars <- diff_names <- c()
  i <- 1
  for (j in 1:4) {
    for (k in 1:4) {
      if (j >= k) next
      score_diff <- rowSums( scoring_function(models[[j]], observations) 
                             - scoring_function(models[[k]], observations))
      means[i] <- mean(score_diff)
      vars[i] <- var(score_diff)
      # Adjust name of difference such that model with larger score comes first
      if (means[i] >= 0) {
        diff_names[i] <- paste(model_names[j], model_names[k], sep = "$-$")
      } else {
        diff_names[i] <- paste(model_names[k], model_names[j], sep = "$-$")
      }
      i <- i + 1
    }
  }
  return(list("means" = means,
              "vars" = vars,
              "diff_names" = diff_names))
}


# Print table of mean and variance of score differences and their
# detectable differences
sample_size_table <- function(diff_means, diff_vars, diff_names, ndays,
                              score_name, file_path, scaling = 1) {
  # Input values:
  # diff_means  - Means of score differences
  # diff_vars   - Variances of score differences
  # diff_names  - LateX strings of score differences
  # ndays       - Number of days in testing period
  # score_name  - Name of the used scoring function
  # file_path   - File path for .tex file
  # scaling     - Factor to scale the table values (optional)
  diff_means <- abs(scaling^2 * diff_means)
  diff_vars <- scaling^2 * diff_vars
  c_string <- paste(rep("c", length(diff_means)), collapse = "")
  begin <- paste0("\\begin{tabular}{l ", c_string, "}")
  mean_order <- order(diff_means, decreasing = TRUE)
  head <- paste(diff_names[mean_order], collapse = " & ")
  head <- paste0(score_name, " & ", head, " \\\\")
  # Write teX code to file
  write(begin, file_path)
  write("\\hline ", file_path, append = T)
  write(head, file_path, append = T)
  write("\\hline ", file_path, append = T)
  # Row with means
  vals <- paste(sprintf("%.3f", diff_means[mean_order]), collapse = " & ")
  vals <- paste0("Mean $m$ & ", vals, " \\\\")
  write(vals, file_path, append = T)
  # Row with variances
  vals <- paste(sprintf("%.3f", diff_vars[mean_order]), collapse = " & ")
  vals <- paste0("Variance $s^2$ &", vals, " \\\\")
  write(vals, file_path, append = T)
  write("\\hline", file_path, append = T)
  # Rows with detectable differences for full and weekly sample
  for (n in c(ndays, floor(ndays/7))) {
    l_est <- lehr_estimate(diff_vars[mean_order], n)
    vals <- paste(sprintf("%.3f", scaling * l_est), collapse = " & ")
    vals <- paste0("$d_{", n, "}$ & ", vals, " \\\\")
    write(vals, file_path, append = T)
  }
  write("\\hline", file_path, append = T)
  write("\\end{tabular}", file_path, append = T)
}


# Plot ACF for score differences of all model combinations
plot_score_diffs_acf <- function(models, observations, scoring_function,
                                 file_path) {
  # Input values:
  # models            - List of model outputs matrices
  # observations      - Observation matrix
  # scoring_function  - Scoring function to compute score differences
  # file_path         - File path for .pdf file
  diff_means <- diff_names <- c()
  diffs <- list()
  i <- 1
  for (j in 1:4) {
    for (k in 1:4) {
      if (j >= k) next
      score_diff <- rowSums( scoring_function(models[[j]], observations) 
                             - scoring_function(models[[k]], observations))
      diffs[[i]] <- score_diff
      diff_means[i] <- mean(score_diff)
      # Adjust name of difference such that model with larger score comes first
      if (diff_means[i] >= 0) {
        diff_names[i] <- paste(model_names[j], model_names[k], sep = "-")
      } else {
        diff_names[i] <- paste(model_names[k], model_names[j], sep = "-")
      }
      i <- i + 1
    }
  }
  # Sort by size of mean score difference (same order as in function
  # sample_size_table for consistency)
  mean_order <- order(abs(diff_means), decreasing = TRUE)
  pdf(file_path, width = 8, height = 5)
  par(mfrow=c(2,3), mar = c(2, 2, 2, 1))
  for (i in mean_order) {
    acf(diffs[[i]], main = "", lag.max = 16, xlab = "", ylab = "")
    title(diff_names[i], line = 0.3)
  }
  dev.off()
}
