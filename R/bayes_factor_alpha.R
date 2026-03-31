#' Bayes Factor for a Point Null on Alpha
#'
#' Computes a Savage-Dickey style Bayes factor comparing a point null
#' \eqn{alpha = alpha0} against a Uniform(0, 1) alternative for the `al`
#' parameter returned by [bayesian_mcmc()].
#'
#' The posterior density at `alpha0` can be estimated either with a kernel
#' density estimate or by a simple local bin count.
#'
#' @param fit Output from [bayesian_mcmc()].
#' @param alpha0 Point null value to evaluate. Must lie in `[0, 1]`.
#' @param method Density estimator to use. One of `"kde"` or `"bin"`.
#' @param eps_bin Half-width of the bin around `alpha0` when
#'   `method = "bin"`.
#'
#' @return A list with the evaluated point `alpha0`, the posterior density
#'   estimate at that point, the Uniform(0,1) prior density, the Bayes factor,
#'   its logarithm, and the posterior samples used in the calculation.
#'
#' @details
#' In the package output, the `al` parameter is stored after the lineage-level
#' `mu_i` parameters and `sig_xi`, so it is extracted from column `N + 2` of
#' `fit$samples`, where `N = fit$N`.
#'
#' The returned `BF12` is computed as posterior density divided by prior
#' density at `alpha0`. Under the Uniform(0,1) prior, the prior density equals
#' 1 everywhere in `[0,1]`.
#'
#' Values larger than 1 indicate more posterior mass near `alpha0` than would
#' be expected under the Uniform(0,1) prior. Values smaller than 1 indicate
#' evidence against the point null.
#'
#' The `"kde"` method is smoother and usually preferable for moderate or large
#' posterior samples. The `"bin"` method is a simpler local estimate that can be
#' useful as a quick sensitivity check.
#'
#' @examples
#' x <- matrix(
#'   c(
#'     1.1, 1.3, 1.5, 0,
#'     0.9, 1.0, 1.2, 1.4,
#'     1.8, 2.0, 0, 0
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' fit <- bayesian_mcmc(
#'   x,
#'   nIter = 200,
#'   burnin = 100,
#'   thin = 10,
#'   verbose = FALSE
#' )
#'
#' bayes_factor_alpha_point_vs_uniform01(fit, alpha0 = 0.5)
#' bayes_factor_alpha_point_vs_uniform01(fit, alpha0 = 0.5, method = "bin")
#' @export
bayes_factor_alpha_point_vs_uniform01 <- function(
  fit,
  alpha0,
  method = c("kde", "bin"),
  eps_bin = 1e-3
) {
  method <- match.arg(method)

  if (!is.list(fit) || is.null(fit$N) || is.null(fit$samples)) {
    stop("`fit` must be a result returned by `bayesian_mcmc()`.")
  }

  if (!is.numeric(alpha0) || length(alpha0) != 1L || is.na(alpha0)) {
    stop("`alpha0` must be a single numeric value.")
  }

  if (!is.numeric(eps_bin) || length(eps_bin) != 1L || is.na(eps_bin) || eps_bin <= 0) {
    stop("`eps_bin` must be a positive numeric scalar.")
  }

  N <- fit$N
  alpha_samples <- fit$samples[, N + 2]

  alpha_samples <- alpha_samples[alpha_samples >= 0 & alpha_samples <= 1]

  if (!length(alpha_samples)) {
    stop("No alpha samples remain in [0,1].")
  }

  if (alpha0 < 0 || alpha0 > 1) {
    stop("`alpha0` must be in [0,1].")
  }

  if (method == "kde") {
    dens <- stats::density(alpha_samples, from = 0, to = 1, n = 2048)
    post_density_at_a0 <- stats::approx(dens$x, dens$y, xout = alpha0, rule = 2)$y
  } else {
    in_bin <- abs(alpha_samples - alpha0) < eps_bin
    post_density_at_a0 <- mean(in_bin) / (2 * eps_bin)
  }

  prior_density_at_a0 <- 1
  BF12 <- post_density_at_a0 / prior_density_at_a0
  logBF12 <- log(BF12)

  list(
    alpha0 = alpha0,
    method = method,
    alphaSamples = alpha_samples,
    posteriorDensityAtAlpha0 = post_density_at_a0,
    priorDensityAtAlpha0 = prior_density_at_a0,
    BF12 = BF12,
    logBF12 = logBF12,
    nSamplesUsed = length(alpha_samples)
  )
}
