#' Run Bayesian MCMC Inference
#'
#' Fits the Bayesian lineage model implemented in C++ and returns posterior
#' samples and summary statistics.
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
#'   excluded from the likelihood.
#' @param alphaMaxOpt Upper bound for the `alpha` parameter during
#'   initialization.
#' @param useMAPinit Reserved argument kept for compatibility with the C++
#'   entry point.
#' @param steps_ Optional scalar or full parameter vector of proposal step
#'   sizes.
#' @param lb_ Optional lower bound vector.
#' @param ub_ Optional upper bound vector.
#' @param par0_ Optional initial parameter vector.
#'
#' @return A list containing posterior draws, log-posterior values, posterior
#'   summaries, acceptance rates, and parameter bounds.
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
