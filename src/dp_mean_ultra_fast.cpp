// Ultra-fast C++ implementation using Nelder-Mead (no gradients needed)
// Optimized for ~100x speedup

#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Fast objective function evaluation with minimal overhead
inline double objective_fast(double mu, double sigma,
                            double a, double b, int n,
                            double scale_mean, double scale_sd,
                            const double* z, double u1, double u2,
                            double mean_obs, double sd_obs) {

  // Generate and truncate data - compute mean inline
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    double x_i = mu + sigma * z[i];
    x_i = (x_i < a) ? a : ((x_i > b) ? b : x_i);  // Fast clamp
    sum += x_i;
  }
  double mean_trunc = sum / n;

  // Compute SD inline
  double sum_sq = 0.0;
  for (int i = 0; i < n; i++) {
    double x_i = mu + sigma * z[i];
    x_i = (x_i < a) ? a : ((x_i > b) ? b : x_i);
    double diff = x_i - mean_trunc;
    sum_sq += diff * diff;
  }
  double sd_trunc = std::sqrt(sum_sq / (n - 1));

  // Add Laplace noise
  double sign_u1 = (u1 > 0.5) ? 1.0 : -1.0;
  double sign_u2 = (u2 > 0.5) ? 1.0 : -1.0;
  double abs_u1 = std::abs(u1 - 0.5);
  double abs_u2 = std::abs(u2 - 0.5);

  double noise_mean = -scale_mean * sign_u1 * std::log(1.0 - 2.0 * abs_u1);
  double noise_sd = -scale_sd * sign_u2 * std::log(1.0 - 2.0 * abs_u2);

  double mean_private = mean_trunc + noise_mean;
  double sd_private = sd_trunc + noise_sd;

  // Squared error
  double dm = mean_obs - mean_private;
  double ds = sd_obs - sd_private;
  return dm * dm + ds * ds;
}

// Ultra-fast Nelder-Mead for 2D box-constrained optimization
inline void nelder_mead_2d(double& x0, double& x1,
                          double l0, double l1, double u0, double u1,
                          double a, double b, int n,
                          double scale_mean, double scale_sd,
                          const double* z, double u1_val, double u2_val,
                          double mean_obs, double sd_obs) {

  // Project onto constraints
  auto project = [&](double& p0, double& p1) {
    p0 = std::max(l0, std::min(u0, p0));
    p1 = std::max(l1, std::min(u1, p1));
  };

  project(x0, x1);

  // Nelder-Mead parameters
  const double alpha = 1.0;    // reflection
  const double gamma = 2.0;    // expansion
  const double rho = 0.5;      // contraction
  const double sigma = 0.5;    // shrink
  const int max_iter = 50;     // Reduced iterations for speed
  const double tol = 1e-6;

  // Initial simplex (3 points for 2D)
  double p0[2] = {x0, x1};
  double p1[2] = {x0 + 0.1, x1};
  double p2[2] = {x0, x1 + 0.1};

  project(p1[0], p1[1]);
  project(p2[0], p2[1]);

  double f0 = objective_fast(p0[0], p0[1], a, b, n, scale_mean, scale_sd,
                             z, u1_val, u2_val, mean_obs, sd_obs);
  double f1 = objective_fast(p1[0], p1[1], a, b, n, scale_mean, scale_sd,
                             z, u1_val, u2_val, mean_obs, sd_obs);
  double f2 = objective_fast(p2[0], p2[1], a, b, n, scale_mean, scale_sd,
                             z, u1_val, u2_val, mean_obs, sd_obs);

  for (int iter = 0; iter < max_iter; iter++) {

    // Sort points
    if (f0 > f1) { std::swap(f0, f1); std::swap(p0[0], p1[0]); std::swap(p0[1], p1[1]); }
    if (f0 > f2) { std::swap(f0, f2); std::swap(p0[0], p2[0]); std::swap(p0[1], p2[1]); }
    if (f1 > f2) { std::swap(f1, f2); std::swap(p1[0], p2[0]); std::swap(p1[1], p2[1]); }

    // Check convergence
    if (f2 - f0 < tol) break;

    // Centroid of best two points
    double c0 = (p0[0] + p1[0]) * 0.5;
    double c1 = (p0[1] + p1[1]) * 0.5;

    // Reflection
    double r0 = c0 + alpha * (c0 - p2[0]);
    double r1 = c1 + alpha * (c1 - p2[1]);
    project(r0, r1);

    double fr = objective_fast(r0, r1, a, b, n, scale_mean, scale_sd,
                               z, u1_val, u2_val, mean_obs, sd_obs);

    if (fr < f1) {
      if (fr < f0) {
        // Expansion
        double e0 = c0 + gamma * (r0 - c0);
        double e1 = c1 + gamma * (r1 - c1);
        project(e0, e1);

        double fe = objective_fast(e0, e1, a, b, n, scale_mean, scale_sd,
                                   z, u1_val, u2_val, mean_obs, sd_obs);

        if (fe < fr) {
          p2[0] = e0; p2[1] = e1; f2 = fe;
        } else {
          p2[0] = r0; p2[1] = r1; f2 = fr;
        }
      } else {
        p2[0] = r0; p2[1] = r1; f2 = fr;
      }
    } else {
      // Contraction
      double ct0, ct1;
      if (fr < f2) {
        ct0 = c0 + rho * (r0 - c0);
        ct1 = c1 + rho * (r1 - c1);
      } else {
        ct0 = c0 + rho * (p2[0] - c0);
        ct1 = c1 + rho * (p2[1] - c1);
      }
      project(ct0, ct1);

      double fct = objective_fast(ct0, ct1, a, b, n, scale_mean, scale_sd,
                                  z, u1_val, u2_val, mean_obs, sd_obs);

      if (fct < f2) {
        p2[0] = ct0; p2[1] = ct1; f2 = fct;
      } else {
        // Shrink
        p1[0] = p0[0] + sigma * (p1[0] - p0[0]);
        p1[1] = p0[1] + sigma * (p1[1] - p0[1]);
        p2[0] = p0[0] + sigma * (p2[0] - p0[0]);
        p2[1] = p0[1] + sigma * (p2[1] - p0[1]);
        project(p1[0], p1[1]);
        project(p2[0], p2[1]);

        f1 = objective_fast(p1[0], p1[1], a, b, n, scale_mean, scale_sd,
                            z, u1_val, u2_val, mean_obs, sd_obs);
        f2 = objective_fast(p2[0], p2[1], a, b, n, scale_mean, scale_sd,
                            z, u1_val, u2_val, mean_obs, sd_obs);
      }
    }
  }

  // Return best point
  x0 = p0[0];
  x1 = p0[1];
}

//' Ultra-fast C++ DP-TOST using Nelder-Mead (100x target)
//'
//' Maximum performance implementation with minimal overhead
// [[Rcpp::export]]
List tost_dp_one_sample_ultra_fast(double a, double b, int n, double epsilon,
                                    double mean_private_obs, double sd_private_obs,
                                    double lower, double upper, int B, double alpha,
                                    int seed) {

  // Set seed
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);

  // Compute scale parameters
  const double delta_f_mean = (b - a) / n;
  const double delta_f_sd = (b - a) / std::sqrt(n - 1);
  const double scale_mean = delta_f_mean / (epsilon * 0.5);
  const double scale_sd = delta_f_sd / (epsilon * 0.5);

  // Generate all random numbers upfront
  NumericVector z_all = Rcpp::rnorm(B * n, 0.0, 1.0);
  NumericVector u_all = Rcpp::runif(B * 2, 0.0, 1.0);

  // Pre-allocate
  NumericVector mu_estimates(B);
  NumericVector sigma_estimates(B);

  // Bounds
  const double lower_mu = a;
  const double lower_sigma = 0.000001;
  const double upper_mu = b;
  const double upper_sigma = (b - a) * 3;

  // Main loop - fully optimized
  for (int i = 0; i < B; i++) {

    // Get pointer to data (zero copy)
    const double* z_i = &z_all[i * n];
    const double u1_i = u_all[i * 2];
    const double u2_i = u_all[i * 2 + 1];

    // Initial guess
    double theta0 = mean_private_obs;
    double theta1 = std::max(lower_sigma, std::min(upper_sigma, sd_private_obs));

    // Optimize
    nelder_mead_2d(theta0, theta1,
                   lower_mu, lower_sigma, upper_mu, upper_sigma,
                   a, b, n, scale_mean, scale_sd,
                   z_i, u1_i, u2_i,
                   mean_private_obs, sd_private_obs);

    mu_estimates[i] = theta0;
    sigma_estimates[i] = theta1;
  }

  // Compute quantiles
  NumericVector mu_sorted = clone(mu_estimates);
  std::sort(mu_sorted.begin(), mu_sorted.end());

  int idx_lower = static_cast<int>(std::floor(alpha * B));
  int idx_upper = static_cast<int>(std::floor((1.0 - alpha) * B));
  idx_lower = std::max(0, std::min(B - 1, idx_lower));
  idx_upper = std::max(0, std::min(B - 1, idx_upper));

  const double ci_lower = mu_sorted[idx_lower];
  const double ci_upper = mu_sorted[idx_upper];

  // TOST decision
  const int decision = (ci_lower > lower && ci_upper < upper) ? 1 : 0;

  return List::create(
    Named("decision") = decision,
    Named("ci_lower") = ci_lower,
    Named("ci_upper") = ci_upper,
    Named("mean_private_obs") = mean_private_obs,
    Named("sd_private_obs") = sd_private_obs,
    Named("mu_estimates") = mu_estimates,
    Named("sigma_estimates") = sigma_estimates
  );
}
