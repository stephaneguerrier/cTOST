// Simple L-BFGS-B implementation for 2D optimization
// Based on the L-BFGS algorithm with box constraints via projection

#ifndef LBFGSB_H
#define LBFGSB_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

class LBFGSB {
private:
    const int max_iter = 200;
    const double tol = 1e-8;
    const int m = 5;  // Memory parameter for L-BFGS

    // Project point onto box constraints
    void project(std::vector<double>& x, const std::vector<double>& lower,
                const std::vector<double>& upper) {
        for (size_t i = 0; i < x.size(); i++) {
            x[i] = std::max(lower[i], std::min(upper[i], x[i]));
        }
    }

    // Backtracking line search
    double lineSearch(std::function<double(const std::vector<double>&)> f,
                     const std::vector<double>& x,
                     const std::vector<double>& grad,
                     const std::vector<double>& dir,
                     const std::vector<double>& lower,
                     const std::vector<double>& upper) {
        const double c1 = 1e-4;
        const double rho = 0.9;
        double alpha = 1.0;

        double f0 = f(x);
        double slope = 0.0;
        for (size_t i = 0; i < x.size(); i++) {
            slope += grad[i] * dir[i];
        }

        for (int iter = 0; iter < 20; iter++) {
            std::vector<double> x_new = x;
            for (size_t i = 0; i < x.size(); i++) {
                x_new[i] += alpha * dir[i];
            }
            project(x_new, lower, upper);

            double f_new = f(x_new);
            if (f_new <= f0 + c1 * alpha * slope) {
                return alpha;
            }
            alpha *= rho;
        }
        return alpha;
    }

public:
    std::vector<double> optimize(
        std::function<double(const std::vector<double>&)> f,
        std::function<void(const std::vector<double>&, std::vector<double>&)> grad_f,
        std::vector<double> x0,
        const std::vector<double>& lower,
        const std::vector<double>& upper) {

        int n = x0.size();
        std::vector<double> x = x0;
        project(x, lower, upper);

        std::vector<double> grad(n);
        grad_f(x, grad);

        // Storage for L-BFGS
        std::vector<std::vector<double>> s_history;
        std::vector<std::vector<double>> y_history;
        std::vector<double> rho_history;

        for (int iter = 0; iter < max_iter; iter++) {
            // Check convergence
            double grad_norm = 0.0;
            for (int i = 0; i < n; i++) {
                grad_norm += grad[i] * grad[i];
            }
            grad_norm = std::sqrt(grad_norm);

            if (grad_norm < tol) {
                break;
            }

            // Compute search direction using L-BFGS
            std::vector<double> dir = grad;
            for (int i = 0; i < n; i++) {
                dir[i] = -grad[i];
            }

            int k = s_history.size();
            if (k > 0) {
                std::vector<double> q = grad;
                std::vector<double> alpha_vals(k);

                // First loop
                for (int i = k - 1; i >= 0; i--) {
                    double alpha_i = 0.0;
                    for (int j = 0; j < n; j++) {
                        alpha_i += s_history[i][j] * q[j];
                    }
                    alpha_i *= rho_history[i];
                    alpha_vals[i] = alpha_i;

                    for (int j = 0; j < n; j++) {
                        q[j] -= alpha_i * y_history[i][j];
                    }
                }

                // Scaling
                double gamma = 1.0;
                if (k > 0) {
                    double sy = 0.0, yy = 0.0;
                    for (int j = 0; j < n; j++) {
                        sy += s_history[k-1][j] * y_history[k-1][j];
                        yy += y_history[k-1][j] * y_history[k-1][j];
                    }
                    if (yy > 0) gamma = sy / yy;
                }

                std::vector<double> r = q;
                for (int j = 0; j < n; j++) {
                    r[j] *= gamma;
                }

                // Second loop
                for (int i = 0; i < k; i++) {
                    double beta = 0.0;
                    for (int j = 0; j < n; j++) {
                        beta += y_history[i][j] * r[j];
                    }
                    beta *= rho_history[i];

                    for (int j = 0; j < n; j++) {
                        r[j] += s_history[i][j] * (alpha_vals[i] - beta);
                    }
                }

                dir = r;
                for (int j = 0; j < n; j++) {
                    dir[j] = -r[j];
                }
            }

            // Line search
            double alpha = lineSearch(f, x, grad, dir, lower, upper);

            // Update x
            std::vector<double> x_new = x;
            for (int i = 0; i < n; i++) {
                x_new[i] += alpha * dir[i];
            }
            project(x_new, lower, upper);

            // Compute new gradient
            std::vector<double> grad_new(n);
            grad_f(x_new, grad_new);

            // Update history
            std::vector<double> s(n), y(n);
            for (int i = 0; i < n; i++) {
                s[i] = x_new[i] - x[i];
                y[i] = grad_new[i] - grad[i];
            }

            double sy = 0.0;
            for (int i = 0; i < n; i++) {
                sy += s[i] * y[i];
            }

            if (sy > 1e-10) {
                if ((int)s_history.size() >= m) {
                    s_history.erase(s_history.begin());
                    y_history.erase(y_history.begin());
                    rho_history.erase(rho_history.begin());
                }
                s_history.push_back(s);
                y_history.push_back(y);
                rho_history.push_back(1.0 / sy);
            }

            x = x_new;
            grad = grad_new;
        }

        return x;
    }
};

#endif
