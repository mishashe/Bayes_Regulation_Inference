\name{bayesian_mcmc}
\alias{bayesian_mcmc}
\title{Run Bayesian MCMC Inference}
\description{
Fits the Bayesian lineage model implemented in C++ and returns posterior
samples and summary statistics.
}
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
\item{vM}{Numeric matrix of observations with one lineage per row and one
time point per column.}

\item{nIter}{Total number of MCMC iterations.}

\item{burnin}{Number of burn-in iterations to discard.}

\item{thin}{Thinning interval for retained draws.}

\item{seed}{Random seed passed into R before sampling.}

\item{adaptEvery}{Adaptation interval during burn-in.}

\item{targetAccept}{Target acceptance rate used for proposal adaptation.}

\item{verbose}{Whether to print periodic progress updates.}

\item{printEvery}{Iteration interval for progress output.}

\item{doMAP}{Reserved compatibility argument forwarded to C++.}

\item{zeroIsPadding}{Whether zeros should be treated as padding values and
excluded from the likelihood. Set this to \code{FALSE} if zero is a real
observed value in your data.}

\item{alphaMaxOpt}{Upper bound used when initializing alpha.}

\item{useMAPinit}{Reserved compatibility argument forwarded to C++.}

\item{steps_}{Optional scalar or vector of proposal step sizes.}

\item{lb_}{Optional lower bound vector of length \code{nrow(vM) + 5}.}

\item{ub_}{Optional upper bound vector of length \code{nrow(vM) + 5}.}

\item{par0_}{Optional initial parameter vector of length \code{nrow(vM) + 5}.}
}
\value{
A list containing posterior draws, posterior summaries, acceptance rates,
final proposal scales, and parameter bounds. Important entries include
\code{samples}, \code{logPost}, \code{mean}, \code{median}, \code{std},
\code{names}, \code{acceptRateOverall}, \code{acceptRatePerParam},
\code{finalSteps}, \code{lb}, \code{ub}, and \code{par0}.
}
\details{
The input matrix \code{vM} is interpreted as one lineage per row and one
time point per column. Missing values (\code{NA} and \code{NaN}) are always
ignored. When \code{zeroIsPadding = TRUE}, zeros are also ignored, which is
useful when lineages with different lengths are stored in the same matrix.

The output parameter names are ordered as:
\describe{
  \item{\code{mu_1}, ..., \code{mu_N}}{Lineage-specific mean parameters.}
  \item{\code{sig_xi}}{Process noise scale.}
  \item{\code{al}}{Autoregressive update parameter used internally as
  \code{1 - phi}.}
  \item{\code{mu0}}{Global mean across lineages.}
  \item{\code{sig0}}{Across-lineage standard deviation.}
  \item{\code{sig_meas}}{Measurement noise scale.}
}

\code{samples} is a matrix with one retained MCMC draw per row. Columns
follow the order given in \code{names}.
}
\examples{
x <- matrix(
  c(
    1.1, 1.3, 1.5, 0,
    0.9, 1.0, 1.2, 1.4,
    1.8, 2.0, 0, 0
  ),
  nrow = 3,
  byrow = TRUE
)

fit <- bayesian_mcmc(
  x,
  nIter = 200,
  burnin = 100,
  thin = 10,
  verbose = FALSE
)

fit$names
fit$mean
}
