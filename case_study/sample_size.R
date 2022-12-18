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
      file_path <- file.path(fpath, paste("sample/pacf", j, k, day_names[i],
                                          ".pdf", sep="_"))
      main = paste0("Partial autocorrelation of score differences between ",
                    mnames[j], " and ", mnames[k], " on ", day_names[i], "s")
      plot_pacf(diffs, file_path, main)
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
    diff_names[ind] <- paste0(mnames[k], " - ", mnames[j])
    ind <- ind + 1
  }
}

file_path <- file.path(fpath, "sample/var_comparison.pdf")
pdf(file_path, width = 8, height = 6)
plot(0, 0, col="white", xlim = c(1, 9), ylim = c(0, 15), ylab = "Varianz",
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
