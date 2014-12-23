# Function for running simulations in both R and MATLAB
# Kelsey Grantham <kelsey.grantham@mq.edu.au>
# 
# run:
#      run.sims('simname')
# where
#      'simname' is the name of a subdirectory that contains setup.R

# Set working directory using: setwd('path_to_folder')

# install modified version of MIICD package
#   fixed error and modified code to construct survival curve at common set of points
install.packages('MIICD_1.1.tar.gz', repos=NULL, type='source')

source('functions.R')

################# Main Simulation Function #########

output.banner <- function(sim.group, sim.name, file.dir, N) {
  cat('+---------------------------------------------------------------------\n')
  cat('| Interval censoring simulation script\n')
  cat('| Kelsey Grantham <kelsey.grantham@mq.edu.au>\n')
  cat('|\n')
  cat(paste('| simulation  =',sim.group,'/',sim.name,'\n'))
  cat(paste('| file.dir    =',file.dir,'\n'))
  cat(paste('| N           =',N,'\n'))
  cat(paste('| started at  =', date(), '\n'))
  cat('+---------------------------------------------------------------------\n')
}

# Get the base directory for the simulation program. Defaults to the 
# current R directory as defined by getwd(), but can be set using the
# SIM_BASE_DIR environment variable
get.base.dir <- function() {
  bdir <- Sys.getenv('SIM_BASE_DIR')
  return(ifelse(bdir == '', getwd(), bdir))
}

# within the simulation group <sim.group>, run the
# simulation set associated with <sim.name>
run.sims <- function(sim.group, sim.name) {
  started.at <- proc.time()
  BASE_DIR <- get.base.dir()
  
  # get config for simulation named sim.name
  file.dir = file.path(BASE_DIR, sim.group, sim.name)
  source(file.path(file.dir, 'setup.R'))

  output.banner(sim.group, sim.name, file.dir, N)

  # step 1: generate the data files
  write.data.files(N, file.dir, seed, generate.data, smooth)
  
  # step 2: simulations in MATLAB and collating simulation files
  if (enable.matlab.sims) {
    run.matlab.sims(N, file.dir, seed, matlab.scripts.dir, matlab.prog)
  }
  
  # step 2b: simluations in Octave on server
  if (enable.octave.sims) {
    run.octave.sims(N, file.dir, seed, matlab.scripts.dir)
  }
  
  # step 3: simulations in R
  if (enable.R.sims) {
    run.R.sims(N, file.dir)
  }
  
  # step 4: collate simulation files
  collate.results(N, file.dir)
  
  # announce we're done. :)
  mn <- (proc.time() - started.at)['elapsed'] / 60
  system(sprintf('say simulations completed in %.1f minutes', mn))
  cat(sprintf('Simulations completed in %.1f minutes\n', mn))
}