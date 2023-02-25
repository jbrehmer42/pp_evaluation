############################################
## Simulation study - Auxiliary functions ##

# Scoring function S1 (intensity)
# Defined in Example 3.5
S1_const <- 0.1
ScoreIntensity1 <- function(fun, norm, ppat) {
  # Input values:
  # fun  - Intensity function
  # norm - Integral of fun
  # ppat - Point pattern (class ppp)
  score <- S1_const * (ppat$n - norm)^2
  if (ppat$n > 0) {
    sspa <- (-1) * log( fun(ppat$x, ppat$y) ) + log(norm)
    score <- score + sum(sspa)
  }
  return(score)
}

# Scoring function S2 (intensity)
# Defined in Proposition 5.1
ScoreIntensity2 <- function(fun, norm, ppat) {
  # Input values:
  # fun  - Intensity function
  # norm - Integral of fun
  # ppat - Point pattern (class ppp)
  score <- norm
  if (ppat$n > 0) {
    sspa <- (-1) * log( fun(ppat$x, ppat$y) )
    score <- score + sum(sspa)
  }
  return(score)
}

# Scoring function for the product density
# Defined in Example B.6
S_const <- 1e-5
ScoreProductDensity <- function(fun, norm, ppat) {
  # Input values:
  # fun  - Function on [0,\infty), called \rho_0 in
  #        Example B.6
  # norm - Integral of the product density which
  #        corresponds to fun
  # ppat - Point pattern (class ppp)
  if (ppat$n > 1) {
    dist <- pairdist(ppat)   # compute all distances
    dist <- as.vector(dist[upper.tri(dist)])
    sspa <- - log( fun(dist) ) + log(norm)
    score <- S_const * ( ppat$n * (ppat$n - 1) - norm )^2 + sum(sspa)
  } else {
    score <- S_const * norm^2
  }
  return(score)
}

# Scoring function for the grid cell forecasts (S_cell)
# Defined in Equation (14)
ScoreCells <- function(cells, counts) {
  # Input values:
  # cells  - List of forecast values for cells
  # counts - List of cell count
  # Each list element contains forecasts and counts
  # for one value of n, i.e. one partition
  ScorePoisson <- function(x,y) sum(x - y * log(x))
  score <- mapply(ScorePoisson, cells, counts)
  return(score)
}


# Compute integral of an intensity function on [0,1]^2
integrate2D <- function(fun, ...) {
  # Input value:
  # fun - Function on [0,1]^2 e.g. intensity function
  f <- function(x) fun(x[1],x[2], ...)
  num <- hcubature(f, c(0,0), c(1,1))
  return(num$integral)
}

# Compute integral of a product density on [0,1]^4 
integrate4D <- function(fun, tol = 1e-5, ...) {
  # Input value:
  # fun - Function on [0,\infty) which defines a product
  #       density on [0,1]^4. This function is called
  #       \rho_0 in Example B.6
  f <- function(x) fun(sqrt( (x[1] - x[3])^2 + (x[2] - x[4])^2 ), ...)
  num <- hcubature(f, c(0,0,0,0), c(1,1,1,1), tol = tol)
  return(num$integral)
}


# Simulate realizations of the point processes
simu <- function(name, N) {
  # This function wraps the simulation routines for the
  # four inhomogeneous and  the three homogeneous point
  # process models
  # Input values:
  # name - Name of point process model for simulation
  # N    - Sample size
  if (name == "Poisson") dat <- rpoispp(f0, nsim = N)
  if (name == "DPP") {
    dat <- simulate( dppGauss(lambda = DPPvar, alpha = DPPscale, d = 2), nsim = N)
    # Thin the homogeneous DPP
    f0thin <- function(ppat) rthin(ppat, function(x,y) f0(x,y)/DPPvar )
    dat <- lapply(dat, f0thin)
  }
  if (name == "LGCP") dat <- rLGCP(LGCPcov, LGCPmu,
                                   param = list(var = LGCPvar, scale = LGCPscale), nsim = N)
  if (name == "Thomas") dat <- rThomas(TPPkappa, TPPscale, TPPmu, nsim = N)
  if (name == "homLGCP") dat <- rLGCP(LGCPcov, LGCPmu,
                                      param = list(var = LGCPvar, scale = LGCPscale), nsim = N)
  if (name == "homPoisson") dat <- rpoispp(lambda, nsim = N)
  if (name == "homDPP") dat <- simulate( dppGauss(lambda = lambda, alpha = DPPscale, d = 2), nsim = N)
  return(dat)
}

# Compute the results of Diebold-Mariano (DM) tests
DMmatrix <- function(s, alpha) {
  # Input values:
  # s     - Array of realized scores. Columns are forecasts,
  #         rows are realizations (nrow is the sample size)
  # alpha - Level of the DM test
  n <- dim(s)[1]
  ncol <- dim(s)[2]
  diffs <- apply(s, 1, function(x) outer(x, x, "-"))
  meandiffs <- rowMeans(diffs)
  vardiffs <- apply(diffs, 1, function(x) var(x))
  # Check whether diffs are significantly far from zero
  logdiffs <- (meandiffs < qnorm(alpha) * sqrt(vardiffs/n))
  # Convert to matrix
  return(matrix(logdiffs, ncol = ncol))
}

# Produces one simulation experiment as described in Section 4
# and Section C of the Supplement
experiment <- function(phi, scoringfun, forecast_functions, forecast_norms, M, N) {
  # Input values:
  # phi        - Name of point process model for simulation
  # scoringfun - Scoring function for intensity measure or product
  #              density
  # forecast_functions - List of forecasts
  # forecast_norms - List of integrals corresponding to forecast_functions
  # M          - Number of repetitions
  # N          - Sample size
  n_forecast <- length(forecast_functions)
  scores <- matrix(0, nrow = M, ncol = n_forecast)
  scores_DM <- array(0, dim = c(n_forecast, n_forecast, M))
  for (j in 1:M) {
    # Simulate N realizations of Phi
    dat <- simu(phi, N)
    s <- matrix(0, nrow = N, ncol = n_forecast)
    # Scores for this sample
    for (i in 1:n_forecast) s[ ,i] <- sapply(dat, FUN = scoringfun, fun = forecast_functions[[i]], norm = forecast_norms[i])
    # Mean scores
    scores[j, ] <- colMeans(s)
    # Restults of DM tests for this sample
    scores_DM[ , ,j] <- DMmatrix(s, alpha)
  }
  # Mean results of DM tests over all repetitions
  DM <- apply(scores_DM, c(1,2), "mean")
  return(list(Scores = scores, DM = DM))
}

# Produce one simulation experiment for the cell approximation
# as described in Section C.1 of the Supplement
experiment_cells <- function(phi, forecast_cell_integrals, n_partition, forecast_base, M, N) {
  # Input values:
  # phi        - Name of point process model for simulation
  # forecast_cell_integrals - list of arrays which contain the cell forecasts
  # n_partition       - number of intervals in the partition of each axis
  # forecast_base      - Index of the forecast for which the convergence of
  #              the preferring probabilities should be studied
  # M          - Number of repetitions
  # N          - Sample size
  n_forecast <- length(forecast_cell_integrals)
  ngrid <- length(n_partition)
  cells_DM <- array(0, dim = c(n_forecast-1, ngrid, M))
  for (j in 1:M) {
    # Simulate N realizations from Phi
    dat <- simu(phi, N)
    s <- array(0, dim = c(N, n_forecast, ngrid))
    nvals <- rep(n_partition, each = length(dat))
    # Count number of points in each grid cell
    count_list <- mapply(cellCount, dat, nvals, USE.NAMES = F)
    count_list <- split(count_list, rep(1:length(dat), ngrid))
    # Scores for this sample
    for (i in 1:n_forecast) {
        s[ , i, ] <- t( sapply(count_list, FUN = ScoreCells, cells = forecast_cell_integrals[[i]]) )
    }
    # Results of DM tests for this sample. Only the DM test which
    # compares the forecast forecast_base to its competitors is saved
    cells_DM[ , ,j] <- apply(s, 3, FUN = function(x) DMmatrix(x, alpha = alpha)[forecast_base, -forecast_base])
  }
  # Mean results of DM tests over all repetitions
  DM <- apply(cells_DM, c(1,2), "mean")
  return(DM)
}


# Print DM table to a .tex file
array2teX <- function(tab, forecast_index, file_path, color = NULL) {
  # Take a matrix and print its values in a tabular
  # environment to a .tex file
  # Input values:
  # tab      - Square matrix
  # forecast_index    - Indices of forecasts for which the DM tests
  #            were performed
  # file_path - File path for .tex file
  # color    - Background color for cells (optional)
  k <- dim(tab)[1]
  tab <- round(tab, digits = 2)
  # header
  cstr <- paste(rep("c", k), collapse = "")
  begin <- paste("\\begin{tabular}{c", cstr, "}", collapse = "")
  forecast_names <- paste0("$f_", forecast_index, "$")
  ftex <- paste0( forecast_names, collapse = " & ")
  ftex <- paste(" & ", ftex, " \\\\", collapse = "")
  # Write to file
  write(begin, file_path)
  write(ftex, file_path, append = T)
  for (i in 1:k) {
    # write values
    row <- sprintf('%.2f', tab[i, ])
    if (is.character(color)) {
      # Add background color according to value in each cell
      row_cols <- paste0("\\cellcolor{", color, "!", round(100 * tab[i, ]), "}" )
      row <- paste0(row_cols, row)
    } 
    row[i] <- " "
    row <- paste(row, collapse = " & ")
    row <- paste0(forecast_names[i], " & ", row, " \\\\")
    write(row, file_path, append = T)
  }
  write("\\end{tabular}", file_path, append = T)
}

ymgp <- c(2,0.6,0)
ycex.axis <- 1
cex.main <- 1.6
cex.axis <- 1.4
gmgp <- c(3, 1.3, 0) 

# Produce boxplots for several experiments
myBoxplot <- function(boxlist, forecast_names, model_names, ylim, file_path, zline = TRUE) {
  # Input values:
  # boxlist  - List of values
  # forecast_names   - Group labels, can be indices of forecasts
  # model_names   - Titles of the individual boxplots
  # ylim     - List of ranges for the boxplots
  # file_path - File path for .pdf file
  # zline    - Plot horizontal line at zero?
  nplots <- length(boxlist)
  n_forecast <- dim(boxlist[[1]])[2]
  nrepeats <- dim(boxlist[[1]])[1]
  if (is.numeric(forecast_names)) forecast_names <- rep(list(paste0("f", forecast_names)), nplots)
  pdf(file_path, width = 14, height = 6)
  par(mfrow = c(1,nplots), cex = 1, cex.axis = cex.axis, cex.main = cex.main,
      mar = c(3.1, 1.6, 2, 0.3), mgp = gmgp)
  for (i in 1:nplots) {
    boxplot(split(boxlist[[i]], rep(1:n_forecast, each = nrepeats)), main = model_names[[i]],
            ylim = ylim[[i]], yaxt = "n", xaxt = "s", col = "lightgray", names = forecast_names[[i]])
    axis(2, mgp = ymgp, cex.axis = ycex.axis)
    if (zline) abline(h = 0)
  }
  dev.off()
}

# Produce 3x2 plot of the intensity function forecasts
plotIntensities <- function(forecast_functions, file_path) {
  # Input values:
  # forecast_functions - List of forecasts
  # file_path   - File path for .pdf file
  # Define grid
  ngrid <- 100
  x <- y <- seq(0, 1, length = ngrid)
  # Compute function values on the grid
  zvals <- sapply(forecast_functions, function(fun) outer(x,y, fun), simplify = "array")
  zlim <- c(min(zvals), max(zvals))
  pdf(file_path, width = 7, height = 10)
  par(mfrow = c(3,2), cex = 1, cex.main = 1, mar = c(2.1, 2, 1.5, 1.5),
      mgp = c(3, 0.7, 0.2), cex.axis = 0.9)
  for (i in 1:6) {
    ftext <- paste("f", i-1, sep = "")
    # Colors
    image(x, y, zvals[ , , i], zlim = zlim, main = paste("Forecast", ftext),
          col = hcl.colors(14, "heat", rev = T), xlab = "", ylab = "", bty = "n")
    # Contour lines
    contour(x, y, zvals[ , , i], col = "gray25", add = TRUE, method = "edge",
            vfont = c("sans serif", "bold"), labcex = 0.8)
  }
  dev.off()
}

# Produce plot of the functions corresponding to the product density 
# forecasts
plotProduct <- function(forecast_functions, file_path) {
  # Input values:
  # forecast_functions - List of forecasts
  # file_path   - File path for .pdf file
  cols <- c("magenta", "green", "black", "blue", "red")
  lwds <- c(3,4,3,3,3)
  ltys <- c("dashed", "dotted", "solid", "dotdash", "A4343434")
  r <- seq(0, 0.2, by = 0.001)
  ylim <- c(0, exp(LGCPvar) * lambda^2)
  pdf(file_path, width = 10, height=6)
  par(cex.axis = cex.axis, cex.main = cex.main,
      mar = c(3, 2.1, 1.5, 0.6) )
  plot(r, forecast_functions[[1]](r), ty="l", col = cols[1], ylim=ylim, lty = ltys[1],
       lwd = lwds[[1]], main = NULL, ylab = NULL)
  for (i in 2:length(forecast_functions)) {
    lines(r, forecast_functions[[i]](r), col = cols[i], lwd = lwds[[i]], lty = ltys[i])
  } 
  legend(x=0.165, y=3150, legend = c("f1", "f2", "f3", "f4", "f5"),
         lwd = lwds, cex = cex.main, lty = ltys, col = cols, seg.len = 3)
  dev.off()
}


# Compute collection of grid cell integrals
cellIntegrals <- function(fun, k) {
  # Input values:
  # fun - Intensity function
  # k   - Number of intervals in the partition of
  #       each axis
  f <- function(x) fun(x[1], x[2])
  lim <- seq(0, 1, len = k+1)
  integrals <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {   # x-axis
    for (j in 1:k) {   # y-axis
      # Integral of fun over grid cell i,j
      num  <- hcubature(f, c(lim[i],lim[j]), c(lim[i+1],lim[j+1]) )
      integrals[i,j] <- num$integral
    }
  }
  return(integrals)
}

# Count the number of points in the grid cells
cellCount <- function(ppat, k) {
  # Intput values:
  # ppat - Point pattern (class ppp)
  # k    - Number of intervals in the partition of
  #        each axis
  lim <- seq(0, 1, len = k+1)
  fleq <- function(a,b) (a <  b)
  fgeq <- function(a,b) (a >= b)
  xleq <- outer(ppat$x, lim[2:(k+1)], fleq)
  xgeq <- outer(ppat$x, lim[1:k], fgeq)
  yleq <- outer(ppat$y, lim[2:(k+1)], fleq)
  ygeq <- outer(ppat$y, lim[1:k], fgeq)
  count <- t( (xleq & xgeq) ) %*% (yleq & ygeq)
  return(count)
}


# Plot the convergence of DM tests preferring frequencies
plotDMconv <- function(DM_cells, DMlims, n_partition, forecast_base, model_names, file_path) {
  # Input values:
  # DM_cells - List of mean DM test results based on S_cell and
  #              different partitions
  # DMlims     - List of mean DM test results based on the scoring
  #              function S2 (function ScoreIntensity2)
  # n_partition       - number of intervals in the partition of each axis
  # forecast_base      - Index of the forecast for which the convergence of
  #              the preferring probabilities should be studied
  # model_names     - Titles of the individual plots
  # file_path   - File path for .pdf file
  cols <- c("black", "magenta", "green", "blue", "red")
  ltys <- c("dashed", "dotted", "dotdash", "twodash", "A4343434")
  lwds <- c(3,4,3,3,3)
  forecast_names <- paste0("f", forecast_index[-forecast_base])
  n_forecast <- length(forecast_index)
  ngrid <- length(n_partition)
  pdf(file_path, width = 14, height=13)
  par(mfrow = c(2,2), cex = 1, cex.axis = cex.axis, cex.main = cex.main,
      mar = c(2.5, 2.5, 3, 0.5), mgp = gmgp)
  for (j in 1:4) {
    plot(1:ngrid, rep(1,ngrid), col = "white", main = model_names[j],
         ylim = c(0,1), xaxt = "n")
    axis(1, 1:ngrid, n_partition)
    lims <- DMlims[[j]][forecast_base, -forecast_base]
    for (i in 1:(n_forecast-1)) {
      lines(1:ngrid, DM_cells[[j]][i, ], col = cols[i], lty = ltys[i], lwd = lwds[i])
      abline(h = lims[i], col = cols[i])
    }
    # Add legend only to first plot for lack of space
    if (j == 1) legend(ngrid - 1.3, 0.38, legend = forecast_names, cex = cex.main,
                       lwd = lwds, col = cols, lty = ltys, seg.len = 3)
  }
  dev.off()
}
