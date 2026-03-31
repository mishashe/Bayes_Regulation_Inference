# BayesRegulationInference

`BayesRegulationInference` is an R package for running Bayesian MCMC inference on lineage-level time series using a compiled C++ backend.

The package exposes one main function, `bayesian_mcmc()`, which takes a numeric matrix of observations and returns posterior samples together with summary statistics and acceptance diagnostics.

## Installation from GitHub

You can install the package directly from your GitHub repository with `remotes`:

```r
install.packages("remotes")
remotes::install_github("mishashe/Bayes_Regulation_Inference")
```

If you prefer `devtools`, this works as well:

```r
install.packages("devtools")
devtools::install_github("mishashe/Bayes_Regulation_Inference")
```

Because the package compiles C++ code, you need a working R build toolchain:

- On Linux: `g++`, `make`, and R development tools.
- On macOS: Xcode Command Line Tools.
- On Windows: Rtools matching your R version.

## Load the package

```r
library(BayesRegulationInference)
```

## Data format

`bayesian_mcmc()` expects a numeric matrix:

- Rows are lineages.
- Columns are time points.
- `NA` and `NaN` values are ignored.
- By default, zeros are treated as padding and ignored.

That default is useful when different lineages have different observed lengths and shorter lineages are padded with zeros to fit into one rectangular matrix.

If zero is a real measurement in your data, set `zeroIsPadding = FALSE`.

## Bundled example input file

The package includes the CSV file `cell_size_simulation.csv` under `extdata`.

You can load it after installation with:

```r
csv_path <- system.file(
  "extdata",
  "cell_size_simulation.csv",
  package = "BayesRegulationInference"
)

vM <- as.matrix(read.csv(csv_path, header = FALSE))
```

## Minimal example

```r
library(BayesRegulationInference)

vM <- matrix(
  c(
    1.10, 1.30, 1.50, 0.00,
    0.95, 1.05, 1.20, 1.35,
    1.80, 2.00, 0.00, 0.00
  ),
  nrow = 3,
  byrow = TRUE
)

fit <- bayesian_mcmc(
  vM,
  nIter = 2000,
  burnin = 1000,
  thin = 10,
  seed = 123,
  verbose = FALSE
)
```

## Inspect the results

The returned object is a list. Common components:

- `samples`: posterior draws matrix.
- `logPost`: log-posterior value for each retained draw.
- `mean`: posterior mean for each parameter.
- `median`: posterior median for each parameter.
- `std`: posterior standard deviation for each parameter.
- `names`: parameter names matching the columns of `samples`.
- `acceptRateOverall`: overall Metropolis acceptance rate.
- `acceptRatePerParam`: per-parameter acceptance diagnostics.
- `finalSteps`: final proposal step sizes.
- `lb`, `ub`, `par0`: bounds and initial values used by the sampler.

Example:

```r
fit$names
fit$mean
fit$std
fit$acceptRateOverall

head(fit$samples)
```

Using the bundled CSV:

```r
csv_path <- system.file(
  "extdata",
  "cell_size_simulation.csv",
  package = "BayesRegulationInference"
)

vM <- as.matrix(read.csv(csv_path, header = FALSE))

fit <- bayesian_mcmc(
  vM,
  nIter = 2000,
  burnin = 1000,
  thin = 10,
  verbose = FALSE
)
```

## Parameter layout

The parameter vector is ordered as:

- `mu_1`, `mu_2`, ..., `mu_N`: lineage-specific mean parameters.
- `sig_xi`: process noise scale.
- `al`: autoregressive update parameter.
- `mu0`: global mean across lineages.
- `sig0`: across-lineage standard deviation.
- `sig_meas`: measurement noise scale.

## Customizing the sampler

You can override defaults when needed:

```r
fit <- bayesian_mcmc(
  vM,
  nIter = 10000,
  burnin = 5000,
  thin = 20,
  adaptEvery = 100,
  targetAccept = 0.30,
  zeroIsPadding = TRUE,
  verbose = TRUE,
  printEvery = 1000
)
```

Advanced arguments:

- `steps_`: scalar or full parameter vector of proposal step sizes.
- `lb_`: custom lower bounds.
- `ub_`: custom upper bounds.
- `par0_`: custom initial parameter values.

The expected length of `lb_`, `ub_`, and `par0_` is `nrow(vM) + 5`.

## Installing from a local clone

If you have already cloned the repository:

```r
install.packages("remotes")
remotes::install_local("Bayes_Regulation_Inference")
```

Or from a shell:

```sh
R CMD INSTALL Bayes_Regulation_Inference
```

## Help inside R

After installation:

```r
?bayesian_mcmc
```
