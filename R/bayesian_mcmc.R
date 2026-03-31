#' Run Bayesian MCMC Inference
#'
#' Fits the Bayesian lineage model implemented in C++ and returns posterior
#' samples and summary statistics.
#'
#' The input matrix `vM` is interpreted as one lineage per row and one time
#' point per column. By default, zero values are treated as padding and ignored,
#' which is useful when lineages have different observed lengths but are stored
#' in a rectangular matrix. Missing values (`NA` and `NaN`) are always ignored.
#'
#' The sampler estimates one lineage-specific mean parameter per row of `vM`,
#' together with global parameters controlling process noise, measurement noise,
#' and shrinkage across lineages.
#'
#' @param vM Numeric matrix of observations. Rows correspond to lineages and
#'   columns correspond to time points.
#' @param nIter Total number of MCMC iterations.
#' @param burnin Number of burn-in iterations to discard.
#' @param thin Thinning interval for retained samples.
#' @param seed Random seed passed into R before sampling.
#' @param adaptEvery Number of burn-in iterations between proposal adaptation
#'   updates.
#' @param targetAccept Target acceptance rate used during adaptation.
#' @param verbose Logical flag controlling progress output.
#' @param printEvery Iteration interval for progress messages when
#'   `verbose = TRUE`.
#' @param doMAP Reserved argument kept for compatibility with the C++ entry
#'   point.
#' @param zeroIsPadding Whether zeros in `vM` should be treated as padding and
#'   excluded from the likelihood. Set this to `FALSE` if zero is a real
#'   observed value in your data.
#' @param alphaMaxOpt Upper bound for the `alpha` parameter during
#'   initialization.
#' @param useMAPinit Reserved argument kept for compatibility with the C++
#'   entry point.
#' @param steps_ Optional scalar or full parameter vector of proposal step
#'   sizes. When a scalar is supplied, the same step size is used for all
#'   parameters.
#' @param lb_ Optional lower bound vector of length `nrow(vM) + 5`.
#' @param ub_ Optional upper bound vector of length `nrow(vM) + 5`.
#' @param par0_ Optional initial parameter vector of length `nrow(vM) + 5`.
#'
#' @return A list containing posterior draws, log-posterior values, posterior
#'   summaries, acceptance rates, and parameter bounds. The main components are:
#'   `samples`, `logPost`, `mean`, `median`, `std`, `names`,
#'   `acceptRateOverall`, `acceptRatePerParam`, `finalSteps`, `lb`, `ub`, and
#'   `par0`.
#'
#' @details
#' Parameter names in the output are ordered as:
#' \describe{
#'   \item{`mu_1`, ..., `mu_N`}{Lineage-specific mean parameters.}
#'   \item{`sig_xi`}{Process noise scale.}
#'   \item{`al`}{Autoregressive update parameter used internally as `1 - phi`.}
#'   \item{`mu0`}{Global mean across lineages.}
#'   \item{`sig0`}{Across-lineage standard deviation.}
#'   \item{`sig_meas`}{Measurement noise scale.}
#' }
#'
#' `samples` is a matrix with one retained MCMC draw per row. Its columns follow
#' the order given in `names`.
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
#' fit$names
#' fit$mean
#' @export
bayesian_mcmc <- function(
  vM,
  nIter = 60000L,
  burnin = 20000L,
  thin = 10L,
  seed = 1L,
  adaptEvery = 50L,
  targetAccept = 0.30,
  verbose = TRUE,
  printEvery = 5000L,
  doMAP = FALSE,
  zeroIsPadding = TRUE,
  alphaMaxOpt = 1.0 - 1e-9,
  useMAPinit = TRUE,
  steps_ = NULL,
  lb_ = NULL,
  ub_ = NULL,
  par0_ = NULL
) {
  if (!is.matrix(vM)) {
    stop("`vM` must be a numeric matrix.")
  }

  storage.mode(vM) <- "double"

  .Call(
    `_BayesRegulationInference__bayesian_mcmc_cpp`,
    vM,
    as.integer(nIter),
    as.integer(burnin),
    as.integer(thin),
    as.integer(seed),
    as.integer(adaptEvery),
    as.numeric(targetAccept),
    as.logical(verbose),
    as.integer(printEvery),
    as.logical(doMAP),
    as.logical(zeroIsPadding),
    as.numeric(alphaMaxOpt),
    as.logical(useMAPinit),
    steps_,
    lb_,
    ub_,
    par0_
  )
}
