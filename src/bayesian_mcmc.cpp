#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>

using namespace Rcpp;

// -------------------------------------------------------------
// Helpers
// -------------------------------------------------------------

static inline bool is_missing_or_nan(double x) {
  return NumericVector::is_na(x) || std::isnan(x);
}

static double vec_mean_std(const std::vector<double>& x, double& s) {
  int n = x.size();
  if (n == 0) {
    s = NA_REAL;
    return NA_REAL;
  }
  double m = std::accumulate(x.begin(), x.end(), 0.0) / n;
  if (n == 1) {
    s = 0.0;
    return m;
  }
  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = x[i] - m;
    ss += d * d;
  }
  s = std::sqrt(ss / (n - 1.0));
  return m;
}

static double vec_mean_only(const std::vector<double>& x) {
  if (x.empty()) return NA_REAL;
  return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

static double vec_sd_only(const std::vector<double>& x) {
  int n = x.size();
  if (n <= 1) return 0.0;
  double m = vec_mean_only(x);
  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = x[i] - m;
    ss += d * d;
  }
  return std::sqrt(ss / (n - 1.0));
}

static double median_stdvec(std::vector<double> x) {
  int n = x.size();
  if (n == 0) return NA_REAL;
  std::sort(x.begin(), x.end());
  if (n % 2 == 1) {
    return x[n / 2];
  } else {
    return 0.5 * (x[n / 2 - 1] + x[n / 2]);
  }
}

static NumericVector col_means_cpp(const NumericMatrix& X) {
  int nrow = X.nrow(), ncol = X.ncol();
  NumericVector out(ncol);
  for (int j = 0; j < ncol; ++j) {
    double s = 0.0;
    for (int i = 0; i < nrow; ++i) s += X(i, j);
    out[j] = s / nrow;
  }
  return out;
}

static NumericVector col_sds_cpp(const NumericMatrix& X, const NumericVector& means) {
  int nrow = X.nrow(), ncol = X.ncol();
  NumericVector out(ncol);
  if (nrow <= 1) {
    std::fill(out.begin(), out.end(), 0.0);
    return out;
  }
  for (int j = 0; j < ncol; ++j) {
    double ss = 0.0;
    for (int i = 0; i < nrow; ++i) {
      double d = X(i, j) - means[j];
      ss += d * d;
    }
    out[j] = std::sqrt(ss / (nrow - 1.0));
  }
  return out;
}

static NumericVector col_medians_cpp(const NumericMatrix& X) {
  int nrow = X.nrow(), ncol = X.ncol();
  NumericVector out(ncol);
  std::vector<double> tmp(nrow);
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) tmp[i] = X(i, j);
    out[j] = median_stdvec(tmp);
  }
  return out;
}

static double objfun_cpp(
    const NumericVector& par,
    int N,
    const std::vector< std::vector<double> >& vCell
) {
  NumericVector mu = par[Range(0, N - 1)];
  double sig_xi      = par[N];
  double alpha       = par[N + 1];
  double mu_epsilon  = par[N + 2];
  double sig_epsilon = par[N + 3];
  double sig_meas    = par[N + 4];

  for (int i = 0; i < par.size(); ++i) {
    if (!R_finite(par[i])) return 1e300;
  }
  if (sig_xi <= 0 || sig_epsilon <= 0 || sig_meas <= 0) return 1e300;
  if (alpha <= 0 || alpha >= 2) return 1e300;

  double phi = 1.0 - alpha;
  double Q   = sig_xi * sig_xi;
  double R   = sig_meas * sig_meas;

  double denom = alpha * (2.0 - alpha);
  if (denom <= 0 || !R_finite(denom)) return 1e300;

  double P0 = Q / denom;
  if (!R_finite(P0) || P0 <= 0) return 1e300;

  double nll = 0.0;
  for (int iLineage = 0; iLineage < N; ++iLineage) {
    const std::vector<double>& y = vCell[iLineage];
    double mu_i = mu[iLineage];

    double a = mu_i;
    double P = P0;

    for (size_t t = 0; t < y.size(); ++t) {
      double S = P + R;
      if (!R_finite(S) || S <= 0) return 1e300;

      double nu = y[t] - a;
      nll += 0.5 * std::log(2.0 * M_PI * S) + 0.5 * (nu * nu) / S;

      double K = P / S;
      double m = a + K * nu;
      double C = (1.0 - K) * P;
      if (C < 1e-15) C = 1e-15;

      a = mu_i + phi * (m - mu_i);
      P = phi * phi * C + Q;
      if (P < 1e-15) P = 1e-15;
    }
  }

  double prior_mu = N * std::log(sig_epsilon);
  double quad = 0.0;
  for (int i = 0; i < N; ++i) {
    double d = mu[i] - mu_epsilon;
    quad += d * d;
  }
  prior_mu += 0.5 * quad / (sig_epsilon * sig_epsilon);

  double nlp = nll + prior_mu;
  if (!R_finite(nlp)) return 1e300;

  return nlp;
}

// [[Rcpp::export]]
Rcpp::List bayesian_mcmc_cpp(
    const NumericMatrix& vM,
    int nIter = 60000,
    int burnin = 20000,
    int thin = 10,
    int seed = 1,
    int adaptEvery = 50,
    double targetAccept = 0.30,
    bool verbose = true,
    int printEvery = 5000,
    bool doMAP = false,
    bool zeroIsPadding = true,
    double alphaMaxOpt = 1.0 - 1e-9,
    bool useMAPinit = true,
    Nullable<NumericVector> steps_ = R_NilValue,
    Nullable<NumericVector> lb_ = R_NilValue,
    Nullable<NumericVector> ub_ = R_NilValue,
    Nullable<NumericVector> par0_ = R_NilValue
) {
  RNGScope scope;
  Function set_seed("set.seed");
  set_seed(seed);

  int N = vM.nrow();
  int T = vM.ncol();

  if (N <= 0 || T <= 0) stop("vM must have positive dimensions.");
  if (nIter <= 0) stop("nIter must be positive.");
  if (thin <= 0) stop("thin must be positive.");
  if (adaptEvery <= 0) stop("adaptEvery must be positive.");

  std::vector< std::vector<double> > vCell(N);
  std::vector<double> vAll;
  vAll.reserve(N * T);

  for (int i = 0; i < N; ++i) {
    std::vector<double> vi;
    vi.reserve(T);

    for (int j = 0; j < T; ++j) {
      double val = vM(i, j);

      if (zeroIsPadding) {
        if (val != 0.0 && !is_missing_or_nan(val)) {
          vi.push_back(val);
          vAll.push_back(val);
        }
      } else {
        if (!is_missing_or_nan(val)) {
          vi.push_back(val);
          vAll.push_back(val);
        }
      }
    }

    if (vi.empty()) stop("Lineage %d has no usable data.", i + 1);
    vCell[i] = vi;
  }

  if (vAll.empty()) stop("No usable observations found in vM.");

  double sAll;
  double mAll = vec_mean_std(vAll, sAll);
  if (!R_finite(sAll) || sAll <= 0) {
    sAll = std::max(1e-6, 0.1 * std::max(1.0, std::fabs(mAll)));
  }

  NumericVector muInit(N);
  for (int i = 0; i < N; ++i) muInit[i] = vec_mean_only(vCell[i]);

  std::vector<double> muInitStd(muInit.begin(), muInit.end());
  double sig0Init;
  double mu0Init = vec_mean_std(muInitStd, sig0Init);
  if (!R_finite(sig0Init) || sig0Init <= 0) {
    sig0Init = std::max(1e-6, 0.5 * sAll);
  }

  double sigXiInit   = std::max(1e-6, 0.5 * sAll);
  double sigMeasInit = std::max(1e-6, 0.5 * sAll);
  double alInit      = 0.5;

  double epsPos = 1e-9;
  double epsA   = 1e-9;

  double dataMin = *std::min_element(vAll.begin(), vAll.end());
  double dataMax = *std::max_element(vAll.begin(), vAll.end());
  double span    = std::max(sAll, (dataMax - dataMin) / 4.0);
  if (!R_finite(span) || span <= 0) {
    span = std::max(1e-3, std::fabs(mAll));
  }

  double lbMu = std::max(epsPos, dataMin - 3.0 * span);
  double ubMu = dataMax + 3.0 * span;
  if (ubMu <= lbMu) ubMu = lbMu + std::max(1.0, std::fabs(lbMu));

  double lbSig = epsPos;
  double ubSig = std::max(10.0 * sAll, 10.0 * epsPos);

  double alphaMax = std::min(alphaMaxOpt, 2.0 - epsA);
  if (alphaMax <= epsA) stop("alphaMax must be > 0 and < 2.");

  int D = N + 5;
  NumericVector lb(D), ub(D), par0(D);

  for (int i = 0; i < N; ++i) {
    lb[i] = lbMu;
    ub[i] = ubMu;
    par0[i] = muInit[i];
  }

  lb[N]     = lbSig;
  lb[N + 1] = epsA;
  lb[N + 2] = lbMu;
  lb[N + 3] = lbSig;
  lb[N + 4] = lbSig;

  ub[N]     = ubSig;
  ub[N + 1] = alphaMax;
  ub[N + 2] = ubMu;
  ub[N + 3] = ubSig;
  ub[N + 4] = ubSig;

  par0[N]     = sigXiInit;
  par0[N + 1] = alInit;
  par0[N + 2] = mu0Init;
  par0[N + 3] = sig0Init;
  par0[N + 4] = sigMeasInit;

  if (lb_.isNotNull()) {
    NumericVector lbIn(lb_);
    if (lbIn.size() != D) stop("lb has wrong length.");
    lb = clone(lbIn);
  }
  if (ub_.isNotNull()) {
    NumericVector ubIn(ub_);
    if (ubIn.size() != D) stop("ub has wrong length.");
    ub = clone(ubIn);
  }
  if (par0_.isNotNull()) {
    NumericVector par0In(par0_);
    if (par0In.size() != D) stop("par0 has wrong length.");
    par0 = clone(par0In);
  }

  for (int j = 0; j < D; ++j) {
    double tiny = 1e-12 * (ub[j] - lb[j]);
    par0[j] = std::min(std::max(par0[j], lb[j] + tiny), ub[j] - tiny);
  }

  int nKeep = (nIter - burnin) / thin;
  if (nKeep <= 0) stop("Need (nIter - burnin) / thin >= 1.");

  NumericVector steps(D);
  for (int j = 0; j < D; ++j) {
    steps[j] = 0.05 * (ub[j] - lb[j]);
    if (steps[j] <= 0) steps[j] = 0.01;
  }

  if (steps_.isNotNull()) {
    NumericVector st(steps_);
    if (st.size() == 1) {
      std::fill(steps.begin(), steps.end(), st[0]);
    } else if (st.size() == D) {
      steps = clone(st);
    } else {
      stop("steps must be scalar or length D.");
    }
  }

  NumericVector cur = clone(par0);
  double curLP = -objfun_cpp(cur, N, vCell);

  NumericVector accCount(D), accTotal(D), tryTotal(D);
  int nAccepted = 0;

  NumericMatrix samples(nKeep, D);
  NumericVector logPost(nKeep);
  int keepIdx = 0;

  for (int it = 1; it <= nIter; ++it) {
    NumericVector prop = clone(cur);

    for (int j = 0; j < D; ++j) {
      double candidate = cur[j] + steps[j] * R::rnorm(0.0, 1.0);
      while (candidate < lb[j] || candidate > ub[j]) {
        candidate = cur[j] + steps[j] * R::rnorm(0.0, 1.0);
      }
      prop[j] = candidate;
    }

    for (int j = 0; j < D; ++j) tryTotal[j] += 1.0;

    double propLP = -objfun_cpp(prop, N, vCell);

    if (std::log(R::runif(0.0, 1.0)) < (propLP - curLP)) {
      cur = prop;
      curLP = propLP;

      for (int j = 0; j < D; ++j) {
        accCount[j] += 1.0;
        accTotal[j] += 1.0;
      }
      nAccepted += 1;
    }

    if (it <= burnin && (it % adaptEvery == 0)) {
      for (int j = 0; j < D; ++j) {
        double ar = accCount[j] / adaptEvery;
        steps[j] = steps[j] * std::exp(ar - targetAccept);
        steps[j] = std::min(std::max(steps[j], 1e-12), 0.5 * (ub[j] - lb[j]));
        accCount[j] = 0.0;
      }
    }

    if (it > burnin && ((it - burnin) % thin == 0)) {
      for (int j = 0; j < D; ++j) samples(keepIdx, j) = cur[j];
      logPost[keepIdx] = curLP;
      keepIdx += 1;
    }

    if (verbose && (it % printEvery == 0)) {
      Rcpp::Rcout << "iter " << it << " / " << nIter
                  << ", kept " << keepIdx << " / " << nKeep
                  << ", acc = " << static_cast<double>(nAccepted) / it
                  << "\n";
    }

    if (it % 10000 == 0) Rcpp::checkUserInterrupt();
  }

  NumericVector postMean   = col_means_cpp(samples);
  NumericVector postMedian = col_medians_cpp(samples);
  NumericVector postStd    = col_sds_cpp(samples, postMean);

  CharacterVector names(D);
  for (int i = 0; i < N; ++i) names[i] = "mu_" + std::to_string(i + 1);
  names[N]     = "sig_xi";
  names[N + 1] = "al";
  names[N + 2] = "mu0";
  names[N + 3] = "sig0";
  names[N + 4] = "sig_meas";

  NumericVector global_mean = postMean[Range(N, N + 4)];
  NumericVector global_std  = postStd[Range(N, N + 4)];

  NumericVector acceptRatePerParam(D);
  for (int j = 0; j < D; ++j) {
    acceptRatePerParam[j] = accTotal[j] / std::max(tryTotal[j], 1.0);
  }
  double acceptRateOverall = static_cast<double>(nAccepted) / nIter;

  return List::create(
    _["samples"] = samples,
    _["logPost"] = logPost,
    _["mean"] = postMean,
    _["median"] = postMedian,
    _["std"] = postStd,
    _["N"] = N,
    _["names"] = names,
    _["global_mean"] = global_mean,
    _["global_std"] = global_std,
    _["acceptRatePerParam"] = acceptRatePerParam,
    _["acceptRateOverall"] = acceptRateOverall,
    _["finalSteps"] = steps,
    _["lb"] = lb,
    _["ub"] = ub,
    _["par0"] = par0
  );
}
