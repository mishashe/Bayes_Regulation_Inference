\name{bayesian_mcmc}
\alias{bayesian_mcmc}
\title{Run Bayesian MCMC Inference}
\usage{
bayesian_mcmc(
  vM,
  nIter = 60000L,
  burnin = 20000L,
  thin = 10L,
  seed = 1L,
  adaptEvery = 50L,
  targetAccept = 0.3,
  verbose = TRUE,
  printEvery = 5000L,
  doMAP = FALSE,
  zeroIsPadding = TRUE,
  alphaMaxOpt = 1 - 1e-09,
  useMAPinit = TRUE,
  steps_ = NULL,
  lb_ = NULL,
  ub_ = NULL,
  par0_ = NULL
)
}
\arguments{
\item{vM}{Numeric matrix of observations with one lineage per row.}

\item{nIter}{Total number of MCMC iterations.}

\item{burnin}{Number of burn-in iterations to discard.}

\item{thin}{Thinning interval for retained draws.}

\item{seed}{Random seed passed into R before sampling.}

\item{adaptEvery}{Adaptation interval during burn-in.}

\item{targetAccept}{Target acceptance rate used for proposal adaptation.}

\item{verbose}{Whether to print periodic progress updates.}

\item{printEvery}{Iteration interval for progress output.}

\item{doMAP}{Reserved compatibility argument forwarded to C++.}

\item{zeroIsPadding}{Whether zeros should be treated as padding values.}

\item{alphaMaxOpt}{Upper bound used when initializing alpha.}

\item{useMAPinit}{Reserved compatibility argument forwarded to C++.}

\item{steps_}{Optional scalar or vector of proposal step sizes.}

\item{lb_}{Optional lower bound vector.}

\item{ub_}{Optional upper bound vector.}

\item{par0_}{Optional initial parameter vector.}
}
\value{
A list containing posterior samples, posterior summaries, acceptance
rates, final proposal scales, and bounds.
}
\description{
Fits the Bayesian lineage model implemented in C++ and returns posterior
samples and summary statistics.
}
