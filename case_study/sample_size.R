###########################
## Sample size calculations

day_names <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
               "Saturday", "Sunday")
day_start <- c(3,4,5,6,7,1,2)

lehr_estimate <- function(var, mean, p) 8 * var / (mean * p)^2

plot_pacf <- function(x, file_path, main) {
  pdf(file_path, width = 8, height = 5.5)
  pacf(x, main = main, lag.max = 7)
  dev.off()
}

plot_acf <- function(x, file_path, main) {
  pdf(file_path, width = 8, height = 5.5)
  acf(x, main = main, lag.max = 7)
  dev.off()
}

plot_diff <- function(x, file_path, main) {
  pdf(file_path, width = 8, height = 5.5)
  plot(1:length(x), x, main = main, xlab = "week", ylab = "score difference")
  dev.off()
}

days_mean <- list()
days_var <- list()

for (i in 1:7) {
  day_ind <- seq(day_start[i], ndays, by = 7)
  # Also consider Squad here
  scores_day <- lapply(models, function(x) rowSums( Spois(x[day_ind, ],
                                                              obs[day_ind, ])))
  day_mean <- day_var <- matrix(0, nrow = nmods, ncol = nmods)
  for (j in 1:4) {
    for (k in 1:4) {
      if (j >= k) next
      diffs <- scores_day[[j]] - scores_day[[k]]
      day_mean[j, k] <- (mean(scores_day[[j]]) + mean(scores_day[[k]]))/2
      day_var[j, k] <- var(diffs)
      file_path <- file.path(fpath, paste("sample/acf", j, k, day_names[i],
                                          ".pdf", sep="_"))
      main = paste0("Autocorrelation of score differences between ",
                    mnames[j], " and ", mnames[k], " on ", day_names[i], "s")
      plot_acf(diffs, file_path, main)
      file_path <- file.path(fpath, paste("sample/diffs", j, k, day_names[i],
                                          ".pdf", sep="_"))
      main = paste0("Score differences between ", mnames[j], " and ", mnames[k],
                    " on ", day_names[i], "s")
      plot_diff(diffs, file_path, main)
    }
  }
  days_mean[[i]] <- day_mean
  days_var[[i]] <- day_var
}

## Plot comparison of variances
days_var2 <- unlist(days_var)
days_var2 <- array(days_var2[days_var2 > 0], dim=c(6, 7))
diff_cols <- c("black", "darkgreen", "blue", "cyan", "red", "magenta")
diff_names <- c()
ind <- 1
for (j in 1:4) {
  for (k in 1:4) {
    if (j >= k) next
    diff_names[ind] <- paste0(mnames[j], " - ", mnames[k])
    ind <- ind + 1
  }
}

file_path <- file.path(fpath, "sample/var_comparison.pdf")
pdf(file_path, width = 8, height = 6)
plot(0, 0, col="white", xlim = c(1, 9), ylim = c(0, 15), ylab = "variance",
     main = "Comparison of score differences variances", xaxt="n", xlab="")
for (i in 1:6) {
  lines(1:7, days_var2[i, ], col = diff_cols[i], lwd = 1.5)
  points(1:7, days_var2[i, ], col = diff_cols[i])
}
par(cex.axis = 0.8)
axis(1, at = 1:7, labels = day_names)
legend("topright", legend = diff_names, col = diff_cols, lty = 1, lwd = 2)
dev.off()


## Plot for every model pair
p <- c(seq(0.02, 0.04, by = 0.001), seq(0.045, 0.1, by = 0.002),
       seq(0.11, 0.3, by = 0.01))

for (j in 1:4) {
  for (k in 1:4) {
    if (j >= k) next
    # Make plot with all days
    file_path <- file.path(fpath, paste("sample/sample_size", mnames[j],
                                        mnames[k], ".pdf", sep = "_"))
    pdf(file_path, width = 8, height =13)
    par(mfrow = c(4, 2), cex = 0.7, mar = c(2, 2, 2, 0.2))
    for (i in 1:7) {
      plot(p, lehr_estimate(days_var[[i]][j, k], days_mean[[i]][j, k], p),
           main = day_names[i], xlab = "relative difference to detect",
           ty = "l", ylab = "number of samples", log = "y")
    }
    dev.off()
  }
}

## Plot for LM-FMC and Friday ([1, 2] and day=5)
p <- c(seq(0.03, 0.04, by = 0.001), seq(0.045, 0.1, by = 0.002),
       seq(0.11, 0.4, by = 0.01))

m1 <- 1
m2 <- 2
day_ind <- 5
file_path <- file.path(fpath, paste0("sample/", "plot_sample_size_pois_",
                                     mnames[m1], "_", mnames[m2], "_",
                                     strtrim(day_names[day_ind], 3), ".pdf"))
pdf(file_path, width = 8, height = 5.5)
par(mar = c(4, 4, 0.5, 0.5))
plot(p, lehr_estimate(days_var[[day_ind]][m1, m2], days_mean[[day_ind]][m1, m2], p),
     main = "", xlab = "relative difference to detect", ty = "l", log = "y",
     ylab = "number of samples", lwd = 2, xaxp = c(0.05, 0.4, 7))
dev.off()


## Print table for supplement
diff_means_pois <- diff_vars_pois <- diff_names_pois <- diff_ess_pois <- c()
diff_means_quad <- diff_vars_quad <- diff_names_quad <- diff_ess_quad <- c()
for (j in 1:4) {
  for (k in 1:4) {
    if (j >= k) next
    score_diff_pois <- rowSums( Spois(models[[j]], obs) - Spois(models[[k]], obs))
    score_diff_quad <- rowSums( Squad(models[[j]], obs) - Squad(models[[k]], obs))
    diff_means_pois <- c(diff_means_pois, mean(score_diff_pois))
    diff_vars_pois <- c(diff_vars_pois, var(score_diff_pois))
    diff_means_quad <- c(diff_means_quad, mean(score_diff_quad))
    diff_vars_quad <- c(diff_vars_quad, var(score_diff_quad))
    if (mean(score_diff_pois) >= 0) {
      diff_names_pois <- c(diff_names_pois, paste(mnames[j], mnames[k], sep = "$-$"))
    } else {
      diff_names_pois <- c(diff_names_pois, paste(mnames[k], mnames[j], sep = "$-$"))
    }
    if (mean(score_diff_quad) >= 0) {
      diff_names_quad <- c(diff_names_quad, paste(mnames[j], mnames[k], sep = "$-$"))
    } else {
      diff_names_quad <- c(diff_names_quad, paste(mnames[k], mnames[j], sep = "$-$"))
    }
    diff_ess_pois <- c(diff_ess_pois, effective_sample_size(score_diff_pois,
                                                            cutoff = 12))
    diff_ess_quad <- c(diff_ess_quad, effective_sample_size(score_diff_quad,
                                                            cutoff = 18))
  }
}

sample_size_table <- function(diff_means, diff_vars, diff_names, score_name,
                              file_path, scaling = 1) {
  diff_means <- abs(diff_means)
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
  n <- 5505
  vals <- paste(sprintf("%.3f", scaling * sqrt(8 * diff_vars[mean_order] / n)),
                collapse = " & ")
  vals <- paste0("$d_{", n, "}$ & ", vals, " \\\\")
  write(vals, file_path, append = T)
  n <- 786
  vals <- paste(sprintf("%.3f", scaling * sqrt(8 * diff_vars[mean_order] / n)),
                collapse = " & ")
  vals <- paste0("$d_{", n, "}$ & ", vals, " \\\\")
  write(vals, file_path, append = T)
  write("\\hline", file_path, append = T)
  write("\\end{tabular}", file_path, append = T)
}

file_path <- file.path(fpath, paste0("sample/", "table_sample_size_pois.tex"))
sample_size_table(diff_means_pois, diff_vars_pois, diff_names_pois,
                  "Poisson score", file_path)
file_path <- file.path(fpath, paste0("sample/", "table_sample_size_quad.tex"))
sample_size_table(100 * diff_means_quad, 100 * diff_vars_quad, diff_names_quad,
                  "Quadratic score", file_path, scaling = 10)


## Plot ACF for daily data
plot_acf_daily <- function(models, obs, scoring_function, file_path) {
  diff_means <- diff_names <- c()
  diffs <- list()
  i <- 1
  for (j in 1:4) {
    for (k in 1:4) {
      if (j >= k) next
      score_diff <- rowSums( scoring_function(models[[j]], obs) 
                             - scoring_function(models[[k]], obs))
      diffs[[i]] <- score_diff
      diff_means[i] <- mean(score_diff)
      if (mean(score_diff) >= 0) {
        diff_names[i] <- paste(mnames[j], mnames[k], sep = "-")
      } else {
        diff_names[i] <-paste(mnames[k], mnames[j], sep = "-")
      }
      i <- i + 1
    }
  }
  mean_order <- order(abs(diff_means), decreasing = TRUE)
  pdf(file_path, width = 8, height = 5)
  par(mfrow=c(2,3), mar = c(2, 2, 2, 1))
  for (i in mean_order) {
    acf(diffs[[i]], main = "", lag.max = 16, xlab = "", ylab = "")
    title(diff_names[i], line = 0.3)
  }
  dev.off()
}

file_path <- file.path(fpath, paste0("sample/", "plot_score_acf_pois", ".pdf"))
plot_acf_daily(models, obs, Spois, file_path)
file_path <- file.path(fpath, paste0("sample/", "plot_score_acf_quad", ".pdf"))
plot_acf_daily(models, obs, Squad, file_path)


## Effective sample size estimation
effective_sample_size <- function(x, cutoff) {
  # Effective sample size (ESS) following Thiebaux and Zwiers (1984)
  rho <- acf(x, lag.max = cutoff, plot = FALSE)$acf
  N <- length(x)
  lags <- 1:cutoff
  denomiator <- 1 + 2 * sum( (1 - lags/N) * rho[lags + 1] )
  return(N / denomiator)
}
