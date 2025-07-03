library(cTOST)

print("1. RUN ALL UNIVARIATE FUNCTIONS")
print("Data: skin")
print(" ")
print(" ")

data(skin)
theta_hat = diff(apply(skin,2,mean))
nu = nrow(skin) - 1
sig_hat = var(apply(skin,1,diff))/nu

# TOST
print("Classical TOST")
print(tost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25)))
print(" ")

# alpha-TOST
print("alpha-TOST")
print(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "alpha"))
print(" ")

# delta-TOST
print("delta-TOST")
print(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "delta"))
print(" ")

# cTOST
print("cTOST (offline)")
print(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "optimal", correction = "offline"))
print(" ")

print("cTOST (no)")
print(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "optimal", correction = "no"))
print(" ")

print("cTOST (bootstrap)")
print(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "optimal", correction = "bootstrap"))
print(" ")

# compare function
print("Compare to TOST:")
compare_to_tost(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "alpha"))
print(" ")
compare_to_tost(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "delta"))
print(" ")
compare_to_tost(ctost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25), method = "optimal", correction = "offline"))
print(" ")

# If the script went this far...
print("All functions ran without error.")



