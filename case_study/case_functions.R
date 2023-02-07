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
load_times <- function(filePath, lastday) {
  # Input values:
  # filePath - File path for the time stamps file
  # lastday  - Last day for which forecasts can be
  #            evaluated.
  # Read time stamps from file
  raw_strings <- read.table(filePath, skip = 1, sep = ";")
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
  matchlastday <- (times$DD == lastday$DD) & (times$MM == lastday$MM) & (times$YY == lastday$YY)
  if (any(matchlastday)) tindex[max(which(matchlastday)+1):length(tindex)] <- FALSE
  times <- times[tindex, ]
  # Return new time stamps and index of corresponding rows
  return(list(times = times, tindex = tindex))
}

# Load forecast values of the models
load_models <- function(filePaths, tindex) {
  # Input values:
  # filePaths - File paths for the model outputs
  # tindex    - Index of model runs (rows) to be used
  # Assume that the file tindex applies to all
  # forecast model files
  # Read forecast values from files
  nmod <- length(filePaths)
  models <- list()
  for (i in 1:nmod) {
    # Rows are days
    # Columns are grid cells
    models[[i]] <- as.matrix( read.table( xzfile(filePaths[i]) ) )
    attr(models[[i]], "dimnames") <- NULL
    # delete rows corresponding to multiple model runs on one day
    models[[i]] <- models[[i]][tindex, ]
  }
  # Return list of model output matrices
  return(models)
}


# Load grid cells (testing region)
load_cells <- function(filePath) {
  # Input values:
  # filePath - File path for the grid cells file
  # Read grid cell file
  cells <- read.csv(filePath, header = FALSE, skip = 1,
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
load_events <- function(filePath, times) {
  # Input values:
  # filePath - File path for the events file
  # times    - Time stamps of testing period
  # Read events file
  events <- read.csv(filePath, header = FALSE, skip = 1,
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
  Tindex <- rep(-1, n)
  for (i in 1:n) {
    DD <- events$DD[i]
    MM <- events$MM[i]
    YY <- events$YY[i]
    gind <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY)
    if (any(gind)) Tindex[i] <- which(gind)
  }
  # Erase events which do not occur during the testing period
  events <- events[(Tindex > 0), ]
  # Add column time index (TI) to the events data frame
  events$TI <- Tindex[Tindex > 0]
  return(events)
}


# Filter out events which do not fall into the 
# testing region. Assign a grid cell number to
# the remaining events
filterRegion <- function(events, cells) {
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
    isLON <- (cell_le < events$LON[i]) & (events$LON[i] <= cell_ri)
    isLAT <- (cell_lo < events$LAT[i]) & (events$LAT[i] <= cell_up)
    if (any(isLON & isLAT) ) {
      # Assign cell number inside testing region
      ind[i] <- cells$N[isLON & isLAT]
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
  obs <- Matrix(0, ncol = ncells, nrow = ndays, sparse = T)
  for (i in 1:ndays) {
    # Collect events in a 7-day period
    ind <- (events$TI >= i) & (events$TI < i + 7)
    if (any(ind)) {
      obs[i, ] <- tabulate(events$N[ind], nbins = ncells)
    } else next
  }
  return(obs)
}

# Print table of mean scores
scores2teX <- function(scores_pois, scores_quad, mnames, filePath) {
  # Input values:
  # scores_pois - Matrix of daily scores
  # scores_quad - Matrix of daily scores
  # mnames      - Names of the forecast models
  # filePath    - File path for the .tex file
  # Take score values and print teX code to
  # file to get values in tabular environment
  # Compute mean scores
  scores_pois <- colMeans(scores_pois)
  scores_quad <- colMeans(scores_quad)
  k <- length(scores_pois)
  begin <- "\\begin{tabular}{l cc}"
  head <- paste("Model", "pois", "quad \\\\", sep = " & ")
  # Write teX code to file
  write(begin, filePath)
  write("\\hline \\hline", filePath, append = T)
  write(head, filePath, append = T)
  write("\\hline", filePath, append = T)
  for (i in 1:k) {
    # Write score values
    row_pois <- sprintf('%.2f', scores_pois[i])
    row_quad <- sprintf('%.4f', scores_quad[i])
    row <- paste(mnames[i], row_pois, row_quad, sep = " & ")
    row <- paste0(row, " \\\\")
    write(row, filePath, append = T)
  }
  write("\\hline", filePath, append = T)
  write("\\end{tabular}", filePath, append = T)
}


# Plot daily scores of the models
plotScores <- function(scores, times, mnames, mcols, filePath, events = NULL, logscale = T) {
  # Input values:
  # scores   - Matrix of daily scores
  # times    - Time stamps of the scores
  # mnames   - Model names for the legend
  # mcols    - Array of colors for the lines
  # filePath - File path for .pdf file
  # events   - Data frame of events (optional)
  # logscale - Should scores be on log scale?
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
  pdf(filePath, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:ndays, rep(0, ndays), ylim = ylim, xlab = "days", col = "white",
       ylab = ylab, main = "", xaxt="n")
  for (i in 1:ncol(scores)) {
    lines(1:ndays, scores[ ,i], col = mcols[i])
  }
  # Add event circles if possible
  if (hasArg(events))  points(evt_days, rep(evt_y, length(evt_days)), cex = 0.8)
  axis(1, at = xats, labels = xlabels)
  legend(-10, ylim[2], mnames, col = mcols, lwd = 2)
  dev.off()
}


# Plot a map of score differences by individually
# coloring the grid cells
diffMap <- function(vals, cells, filePath, events = NULL, borders = T) {
  # Input values:
  # vals     - Numeric values corresponding to
  #            the grid cells
  # cells    - Data frame of grid cells
  # filePath - File path for .pdf file
  # events   - Data frame of events (optional)
  # borders  - Plot national borders? Requires 
  #            R-package "maps"
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
  pdf(filePath, width = 5, height = 4.4)
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
  #         the neighborhodd for aggregation
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

calculate_means_and_vars <- function(models, obs) {
  means_pois <- vars_pois <- names_pois <- c()
  means_quad <- vars_quad <- names_quad <- c()
  i <- 1
  for (j in 1:4) {
    for (k in 1:4) {
      if (j >= k) next
      score_diff_pois <- rowSums( Spois(models[[j]], obs) - Spois(models[[k]], obs))
      score_diff_quad <- rowSums( Squad(models[[j]], obs) - Squad(models[[k]], obs))
      means_pois[i] <- mean(score_diff_pois)
      vars_pois[i] <- var(score_diff_pois)
      means_quad[i] <- mean(score_diff_quad)
      vars_quad[i] <- var(score_diff_quad)
      if (means_pois[i] >= 0) {
        names_pois[i] <- paste(mnames[j], mnames[k], sep = "$-$")
      } else {
        names_pois[i] <- paste(mnames[k], mnames[j], sep = "$-$")
      }
      if (means_quad[i] >= 0) {
        names_quad[i] <- paste(mnames[j], mnames[k], sep = "$-$")
      } else {
        names_quad[i] <- paste(mnames[k], mnames[j], sep = "$-$")
      }
      i <- i + 1
    }
  }
  return(list("means_pois" = means_pois,
              "vars_pois" = vars_pois,
              "names_pois" = names_pois,
              "means_quad" = means_quad,
              "vars_quad" = vars_quad,
              "names_quad" = names_quad))
}

sample_size_table <- function(diff_means, diff_vars, diff_names, ndays,
                              score_name, file_path, scaling = 1) {
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
  vals <- paste(sprintf("%.3f", diff_means[mean_order]), collapse = " & ")
  vals <- paste0("Mean $m$ & ", vals, " \\\\")
  write(vals, file_path, append = T)
  vals <- paste(sprintf("%.3f", diff_vars[mean_order]), collapse = " & ")
  vals <- paste0("Variance $s^2$ &", vals, " \\\\")
  write(vals, file_path, append = T)
  write("\\hline", file_path, append = T)
  n <- ndays
  l_est <- lehr_estimate(diff_vars[mean_order], n)
  vals <- paste(sprintf("%.3f", scaling * l_est), collapse = " & ")
  vals <- paste0("$d_{", n, "}$ & ", vals, " \\\\")
  write(vals, file_path, append = T)
  n <- floor(ndays/7)
  l_est <- lehr_estimate(diff_vars[mean_order], n)
  vals <- paste(sprintf("%.3f", scaling * l_est), collapse = " & ")
  vals <- paste0("$d_{", n, "}$ & ", vals, " \\\\")
  write(vals, file_path, append = T)
  write("\\hline", file_path, append = T)
  write("\\end{tabular}", file_path, append = T)
}
