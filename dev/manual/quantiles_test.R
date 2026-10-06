n_x = 10000
n_y = 20000

nu_x = n_x - 1
nu_y = n_y - 1

sig_x = 0.04
sig_y = 0.04

mu_x = 0.24
mu_y = 0.24

Delta_l = -0.1
Delta_u =  0.1

pi_x = 0.2

alpha = 0.05

delta_l = pnorm(pi_x + Delta_l)
delta_u = pnorm(pi_x + Delta_u)

center = (delta_l + delta_u)/2

(delta = delta_u - center)


x_bar = mu_x + rnorm(1)*sig_x/sqrt(n_x)
y_bar = mu_y + rnorm(1)*sig_y/sqrt(n_y)

sig_hat_x = sqrt(sig_x^2 / nu_x * rchisq(1, nu_x))
sig_hat_y = sqrt(sig_y^2 / nu_y * rchisq(1, nu_y))

theta_hat = (x_bar - y_bar) / sig_hat_y + sig_hat_x / sig_hat_y * pnorm(pi_x)
theta_hat_star = theta_hat - center

l = n_y / n_x
gamma = sig_hat_y^2 / sig_hat_x^2

sigma2_a = 1 + theta_hat^2 / 2 + l/gamma * (1 + pnorm(pi_x)^2/2)

# CI
(CI = theta_hat_star + c(-1, 1) * qnorm(1 - alpha) * sqrt(sigma2_a)/sqrt(n_y))

# Check equivalence
max(abs(CI)) < delta


