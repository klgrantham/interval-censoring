# This file contains functions used to compute performance measures:
#  bias, mean square error and Monte Carlo standard deviation.

# Compute bias
bias <- function(sim.estimates, true.params) {
  bias.sim.estimates <- colMeans(sim.estimates) - true.params
  return(bias.sim.estimates)
}

# Compute MSE
mse <- function(sim.estimates, true.params) {
  sim.estimates.centered <- sweep(as.matrix(sim.estimates), 2, true.params, '-')
  mse.sim.estimates <- apply(sim.estimates.centered^2, 2, mean)
  return(mse.sim.estimates)
}

# Compute Monte Carlo standard deviation
MC.sd <- function(sim.estimates) {
  MC.sd.sim.estimates <- apply(as.matrix(sim.estimates), 2, sd)
  return(MC.sd.sim.estimates)
}
