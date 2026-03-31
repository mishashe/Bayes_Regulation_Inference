\name{bayes_factor_alpha_point_vs_uniform01}
\alias{bayes_factor_alpha_point_vs_uniform01}
\title{Bayes Factor for a Point Null on Alpha}
\usage{
bayes_factor_alpha_point_vs_uniform01(
  fit,
  alpha0,
  method = c("kde", "bin"),
  eps_bin = 0.001
)
}
\arguments{
\item{fit}{Output from \code{\link{bayesian_mcmc}}.}

\item{alpha0}{Point null value to evaluate. Must lie in \code{[0, 1]}.}

\item{method}{Density estimator to use. One of \code{"kde"} or
\code{"bin"}.}

\item{eps_bin}{Half-width of the bin around \code{alpha0} when
\code{method = "bin"}.}
}
\value{
A list containing the evaluated point \code{alpha0}, the posterior density
estimate at that point, the Uniform(0,1) prior density, the Bayes factor,
its logarithm, and the posterior samples used in the calculation.
}
\description{
Computes a Savage-Dickey style Bayes factor comparing a point null
\code{alpha = alpha0} against a Uniform(0,1) alternative for the
\code{al} parameter returned by \code{\link{bayesian_mcmc}}.
}
\details{
The posterior density at \code{alpha0} can be estimated either with a kernel
density estimate or by a simple local bin count.

In the package output, the \code{al} parameter is stored after the
lineage-level \code{mu_i} parameters and \code{sig_xi}, so it is extracted
from column \code{N + 2} of \code{fit$samples}, where \code{N = fit$N}.

The returned \code{BF12} is computed as posterior density divided by prior
density at \code{alpha0}. Under the Uniform(0,1) prior, the prior density
equals 1 everywhere in \code{[0,1]}.

Values larger than 1 indicate more posterior mass near \code{alpha0} than
would be expected under the Uniform(0,1) prior. Values smaller than 1
indicate evidence against the point null.

The \code{"kde"} method is smoother and usually preferable for moderate or
large posterior samples. The \code{"bin"} method is a simpler local estimate
that can be useful as a quick sensitivity check.
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

bayes_factor_alpha_point_vs_uniform01(fit, alpha0 = 0.5)
bayes_factor_alpha_point_vs_uniform01(fit, alpha0 = 0.5, method = "bin")
}
