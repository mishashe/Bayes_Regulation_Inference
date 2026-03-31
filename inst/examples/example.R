library(BayesRegulationInference)

t0 <- proc.time()

csv_path <- system.file(
  "extdata",
  "cell_size_simulation.csv",
  package = "BayesRegulationInference"
)

if (csv_path == "") {
  stop("Could not locate bundled file 'cell_size_simulation.csv'.")
}

vM <- as.matrix(read.csv(csv_path, header = FALSE))

seed <- 1L
nIter <- 10000L
burnin <- as.integer(0.7 * nIter)
thin <- 1L

out <- bayesian_mcmc(
  vM = vM,
  nIter = nIter,
  burnin = burnin,
  thin = thin,
  seed = seed,
  verbose = FALSE
)

N <- out$N
samples <- out$samples
theta <- samples[, (N + 1):(N + 5), drop = FALSE]
param_names <- c("sig_xi", "al", "mu0", "sig0", "sig_meas")
colnames(theta) <- param_names

post_mean <- colMeans(theta)
ci95 <- apply(theta, 2, quantile, probs = c(0.025, 0.975))

cat("\nPosterior summaries (mean +/- half-width of 95% CrI)\n")
cat("---------------------------------------------------\n")
for (j in seq_along(param_names)) {
  lo <- ci95[1, j]
  hi <- ci95[2, j]
  hw <- (hi - lo) / 2
  cat(sprintf(
    "%-8s : %.4f +/- %.4f   (95%% CrI: [%.4f, %.4f])\n",
    param_names[j], post_mean[j], hw, lo, hi
  ))
}

bf_alpha_05 <- bayes_factor_alpha_point_vs_uniform01(
  out,
  alpha0 = 0.5,
  method = "kde"
)

cat("\nBayes factor for al = 0.5 versus Uniform(0,1)\n")
cat("---------------------------------------------\n")
cat(sprintf("BF12   : %.4f\n", bf_alpha_05$BF12))
cat(sprintf("logBF12: %.4f\n", bf_alpha_05$logBF12))

output_dir <- file.path(getwd(), "example_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

summary_df <- data.frame(
  parameter = param_names,
  posterior_mean = as.numeric(post_mean),
  q2.5 = as.numeric(ci95[1, ]),
  q97.5 = as.numeric(ci95[2, ])
)
write.csv(
  summary_df,
  file = file.path(output_dir, "posterior_summary.csv"),
  row.names = FALSE
)

pdf(
  file = file.path(output_dir, "posterior_histograms.pdf"),
  width = 11,
  height = 8.5
)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (j in seq_along(param_names)) {
  hist(
    theta[, j],
    breaks = 40,
    probability = TRUE,
    main = paste("Posterior of", param_names[j]),
    xlab = param_names[j],
    ylab = "Density",
    col = "grey85",
    border = "white"
  )
  abline(v = post_mean[j], col = "blue", lwd = 2)
  abline(v = ci95[, j], col = "red", lty = 2, lwd = 2)
}

plot.new()
legend(
  "center",
  legend = c("Posterior mean", "95% credible interval"),
  col = c("blue", "red"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

par(old_par)
dev.off()

png(
  filename = file.path(output_dir, "alpha_density_kde.png"),
  width = 1200,
  height = 800,
  res = 150
)
dens <- density(theta[, "al"], from = 0, to = 1, n = 2048)
plot(
  dens,
  main = "Posterior density of al",
  xlab = "al",
  ylab = "Density",
  lwd = 2
)
abline(v = 0.5, col = "black", lty = 2, lwd = 2)
abline(v = post_mean["al"], col = "blue", lwd = 2)
legend(
  "topright",
  legend = c("alpha0 = 0.5", "Posterior mean"),
  col = c("black", "blue"),
  lty = c(2, 1),
  lwd = 2,
  bty = "n"
)
dev.off()

cat("\nSaved outputs in:\n")
cat(output_dir, "\n")
cat("- posterior_histograms.pdf\n")
cat("- alpha_density_kde.png\n")
cat("- posterior_summary.csv\n")

print(proc.time() - t0)
