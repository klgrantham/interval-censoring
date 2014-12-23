# setup file common for all of the simulations

# True/false indicator for running MATLAB/Octave and/or R simulations
enable.matlab.sims <- FALSE
enable.octave.sims <- FALSE
enable.R.sims <- TRUE

N <- 10
seed <- 123
beta1 <- 1.5
beta2 <- 0.5
beta3 <- -0.1

# this function gets passed to the data generation function
generate.data <- function() {
  wblcen(samplesize, c(beta1,beta2,beta3), 3, 2, epct)
}

# smoothing parameter
smooth <- 1e-5

# location of MPL scripts
matlab.scripts.dir <- '/Users/kelseygrantham/Documents/MQ_S2_2014/Research/MATLAB'

# path to matlab executable (on Windows, remember to use \\ for /)
matlab.prog <- '/Applications/MATLAB/MATLAB_R2014a.app/bin/matlab'