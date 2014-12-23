# Functions required for running simulations

############### Functions to Generate Data #############

# function to generate data
write.data.files <- function(N, data_dir, myseed, generate.data, smooth) {
  cat(paste('*** Generating',N,'replications to',data_dir,'\n'))
  
  set.seed(myseed)
  
  suppressMessages(library(R.matlab))
  source('wblcen.R')
  
  for (i in 1:N) {
    infile_mat <- file.path(data_dir, sprintf('input_%06d.mat', i));
    infile_rda <- file.path(data_dir, sprintf('input_%06d.rda', i));
    cat(paste('Writing to files',infile_mat,'&',infile_rda,'\n'));
    
    data <- generate.data()
    data$smooth <- smooth # smoothing parameter
    
    writeMat(infile_mat, data=data) # matlab format
    save(data, file=infile_rda)     # R format
  }
  cat('Finished writing files.\n')
}

################ Simulation Functions ############

run.matlab.sims <- function(N, file.dir, myseed, matlab.scripts.dir, matlab.prog) {
  current.directory <- getwd()  
  setwd(matlab.scripts.dir)
  
  # open matlab and execute program to run MPL method
  matlab_params <- paste0(N,",","'",file.dir,"'",",", myseed) # comma-separated list of matlab parameters
  command <- paste0(matlab.prog,' -nodesktop -nosplash -r "run_matlab_sims_mod(',matlab_params,'); quit;"')
  system(command)
  
  setwd(current.directory)
}

run.octave.sims <- function(N, file.dir, myseed, matlab.scripts.dir) {
  # matlab.scripts.dir gets set in common_setup file within simgroup folder
  # change to '/home/kelsey/MATLAB' for running on fisher1 server
  current.directory <- getwd()  
  setwd(matlab.scripts.dir)
  
  matlab_params <- paste0(N,",","'",file.dir,"'",",", myseed) # comma-separated list of matlab parameters
  command <- paste0('octave --quiet --eval "run_matlab_sims_mod(',matlab_params,');"')
  system(command)
  
  setwd(current.directory)
}

# function to store common endpoints across N input data files
common_ymin_ymax <- function(N, data_dir) {
  ymins <- matrix(0,N,1)
  ymaxs <- matrix(0,N,1)
  for (i in 1:N) {
    infile <- file.path(data_dir, sprintf('input_%06d.rda', i));
    load(file=infile)
    rownames(data$y) <- seq_len(nrow(data$y))
    I <-  data$y[,1] != data$y[,2] & data$y[,2] != Inf
    data2 <- data$y[I,]
    t1 <- sort(unique(c(data2[,1] , data2[,2])))
    ymins[i] <- min(t1)
    ymaxs[i] <- max(t1)
  }
  maxymins <- max(ymins)
  minymaxs <- min(ymaxs)
  endpts <- c(maxymins, minymaxs)
}

# function to run R simulations from stored data files
run.R.sims <- function(N, data_dir) {
  cat(paste('*** running ',N,'replications on data in',data_dir,'\n'))
  
  suppressMessages(require(MIICD))
  suppressMessages(require(intcox))
  suppressMessages(require(boot))
  
  endpts <- common_ymin_ymax(N, data_dir)
  
  # Bootstrapping for regression coefficient standard errors for intcox package
  intcox.boot <- function (formula, data, indices) {
    d <- data[indices,]
    fit <- intcox(formula, data=d)
    return(fit$coef)
  }
  
  for (i in 1:N) {
    # INTCOX method
    infile <- file.path(data_dir, sprintf('input_%06d.rda', i));
    outfile <- file.path(data_dir, sprintf('package_intcox_output_%06d.rda', i));
    cat(paste('Processing ',infile,'->',outfile,'...\n'));
    load(file=infile)
    
    data.full <- data.frame(cbind(data$y, data$X))
    data.full.na <- within(data.full, tr <- ifelse(tr==Inf, NA, tr))
    
    cox <- intcox(Surv(tl,tr,type="interval2")~x1+x2+x3, data=data.full.na)
    intcox.betas <- cox$coef
    common_xs <- seq(from=endpts[1], to=endpts[2], length.out=100)
    lambda0mod <- approx(x=cox$time.point, y=cox$lambda0, xout=common_xs, method="linear")
    intcox_bch <- cbind(common_xs, lambda0mod$y)
    
    f <- as.formula(paste("Surv(", paste(colnames(data.full.na)[1], ",", colnames(data.full.na)[2]), ",",
                          "type='interval2')", "~", paste(colnames(data.full.na)[-c(1,2)], collapse=" + ")))
    res <- boot(data = data.full.na, statistic = intcox.boot, R = 10, formula=f) # number of replications set to 10 for quick example
    intcox.se.betas <- apply(res$t, 2, sd)
    
    results <- list(b=cbind(intcox.betas, intcox.se.betas), bch=intcox_bch)
    save(results, file=outfile)
  }
  
  # MIICD method
  for (i in 1:N) {
    infile <- file.path(data_dir, sprintf('input_%06d.rda', i));
    outfile <- file.path(data_dir, sprintf('package_miicd_output_%06d.rda', i));
    cat(paste('Processing ',infile,'->',outfile,'...\n'));
    load(file=infile)
    
    data.full <- data.frame(cbind(data$y, data$X))
    # MIICD.coxph function needs an explanatory variable called "treatment"
    #   which it uses to initialize the estimates by running a crude Cox PL
    #   with "treatment" as the only covariate.
    # We will take the first covariate to be the "treatment" for this purpose
    colnames(data.full)[1:3] <- c("left", "right", "treatment")
    f <- as.formula(paste("~", paste(colnames(data.full)[-c(1,2)], collapse=" + ")))
    sink(file='/dev/null')
    res <- NULL
    try({
      # Note: sometimes fails, in which case res is NULL
      res <- MIICD.coxph(formula = f, data=data.full, method='ANDA', endpts=endpts)
    })
    sink()
    results <- list(b = res$df, surv = res$survmod)
    save(results, file=outfile)
  }
  
  cat('Finished writing files.\n')
}

########### Functions for Collating Results ########

# retrieve results for nth phmpl simulation
get.nth.phmpl.result <- function(n, data_dir) {
  suppressMessages(require(R.matlab))
  phmpl_file <- file.path(data_dir, sprintf('phmpl_ic_ms_output_%06d.mat', n));
  phmpl <- readMat(phmpl_file)
  return(phmpl)
}

# retrieve results for nth MIICD simulation
get.nth.miicd.result <- function(n, data_dir) {
  miicd_file <- file.path(data_dir, sprintf('package_miicd_output_%06d.rda', n));
  load(file=miicd_file)
  list(
    b=results$b,
    surv=results$surv
  )
}

# retrieve results for nth INTCOX simulation
get.nth.intcox.result <- function(n, data_dir) {
  intcox_file <- file.path(data_dir, sprintf('package_intcox_output_%06d.rda', n));
  load(file=intcox_file)
  list(
    b=results$b,
    bch=results$bch
  )
}

# Write results to a table with given file name
write.results <- function(results, data_dir, file.name) {
  write.table(format(results, digits=4), quote=F, file=file.path(data_dir, file.name))
}

# Collate results for beta and hazard curve simulation estimates and compute performance measures
collate.results <- function(N, data_dir) {
  cat(paste('*** tabulating',N,'replications on data in',data_dir,'\n'))
  
  source('performance_measures.R') # loads bias, MSE and Monte Carlo standard deviation functions
  
  # Collate results from MATLAB and compute performance measures
  if (enable.matlab.sims | enable.octave.sims) {
    # data.frame to accumulate the beta estimates and standard deviations
    phmpl_betas <- data.frame()
    phmpl_betas_sd <- data.frame()
    phmpl_bhmod <- data.frame()
    phmpl_bhmod_sd <- data.frame()
    phmpl_bchmod <- data.frame()
    phmpl_bchmod_sd <- data.frame()
    
    for (i in 1:N) {
      cat(ifelse(i %% 50 == 0, '.\n', '.'))
      
      # process ith phmpl replication
      phmpl <- get.nth.phmpl.result(i, data_dir)
      phmpl_betas <- rbind(phmpl_betas, phmpl$b[,1])
      phmpl_betas_sd <- rbind(phmpl_betas_sd, phmpl$b[,2])
      phmpl_bhmod <- rbind(phmpl_bhmod, phmpl$bhmod[,2])
      phmpl_bhmod_sd <- rbind(phmpl_bhmod_sd, phmpl$bhmod[,3])
      phmpl_bchmod <- rbind(phmpl_bchmod, phmpl$bchmod[,1])
      phmpl_bchmod_sd <- rbind(phmpl_bchmod_sd, phmpl$bchmod[,2])
    }
    colnames(phmpl_betas) <- c('beta_1','beta_2','beta_3')
    save(phmpl_betas, file=file.path(data_dir, 'phmpl_betas.rda'))
    colnames(phmpl_betas_sd) <- c('beta_1','beta_2','beta_3')
    save(phmpl_betas_sd, file=file.path(data_dir, 'phmpl_betas_sd.rda'))
    # save baseline hazard values
    phmpl_bhmod_xs <- phmpl$bhmod[,1]
    colnames(phmpl_bhmod) <- seq_len(ncol(phmpl_bhmod))
    phmpl_bh <- list(time=phmpl_bhmod_xs, bh=phmpl_bhmod)
    save(phmpl_bh, file=file.path(data_dir, 'phmpl_bh.rda'))
    colnames(phmpl_bhmod_sd) <- seq_len(ncol(phmpl_bhmod_sd))
    save(phmpl_bh_sd=phmpl_bhmod_sd, file=file.path(data_dir, 'phmpl_bh_sd.rda'))
    colnames(phmpl_bchmod) <- seq_len(ncol(phmpl_bchmod))
    phmpl_bch <- list(time=phmpl_bhmod_xs, bch=phmpl_bchmod)
    save(phmpl_bch, file=file.path(data_dir, 'phmpl_bch.rda'))
    colnames(phmpl_bchmod_sd) <- seq_len(ncol(phmpl_bchmod_sd))
    save(phmpl_bch_sd=phmpl_bchmod_sd, file=file.path(data_dir, 'phmpl_bch_sd.rda'))
    
    # Compute bias
    phmpl_betas_bias <- bias(phmpl_betas, c(beta1, beta2, beta3))
    write.results(results=phmpl_betas_bias, data_dir=data_dir, file.name='phmpl_betas_bias')
    
    # Compute Monte Carlo standard deviation
    phmpl_betas_MC_sd <- MC.sd(phmpl_betas)
    write.results(results=phmpl_betas_MC_sd, data_dir=data_dir, file.name='phmpl_betas_MC_sd')
    
    # Compute average asymptotic standard deviation
    phmpl_betas_mean_sd <- colMeans(phmpl_betas_sd)
    write.results(results=phmpl_betas_mean_sd, data_dir=data_dir, file.name='phmpl_betas_mean_sd')
    
    # Compute mean squared error
    phmpl_betas_mse <- mse(phmpl_betas, c(beta1, beta2, beta3))
    write.results(results=phmpl_betas_mse, data_dir=data_dir, file.name='phmpl_betas_mse')
  }
  
  if (enable.R.sims) {
    # data.frames to accumulate the estimates  
    intcox_betas <- data.frame()
    intcox_betas_sd <- data.frame()
    intcox_bch <- data.frame()
    
    # Collate intcox results
    for (i in 1:N) {
      cat(ifelse(i %% 50 == 0, '.\n', '.'))
      
      # process ith intcox replication
      cox <- get.nth.intcox.result(i, data_dir)
      intcox_betas <- rbind(intcox_betas, cox$b[,1])
      intcox_betas_sd <- rbind(intcox_betas_sd, cox$b[,2])
      intcox_bch <- rbind(intcox_bch, cox$bch[,2])
    }
    colnames(intcox_betas) <- c('beta_1','beta_2','beta_3')
    save(intcox_betas=intcox_betas, file=file.path(data_dir, 'intcox_betas.rda'))
    colnames(intcox_betas_sd) <- c('beta_1','beta_2','beta_3')
    save(intcox_betas_sd=intcox_betas_sd, file=file.path(data_dir, 'intcox_betas_sd.rda'))
    # save baseline cumulative hazard values
    intcox_bch_xs <- cox$bch[,1]
    colnames(intcox_bch) <- seq_len(ncol(intcox_bch))
    intcox_bch_join <- list(time=intcox_bch_xs, cbh=intcox_bch)
    save(intcox_bch_join, file=file.path(data_dir, 'intcox_bch.rda'))
    
    # Compute bias
    intcox_betas_bias <- bias(intcox_betas, c(beta1, beta2, beta3))
    write.results(results=intcox_betas_bias, data_dir=data_dir, file.name='intcox_betas_bias')
    
    # Compute Monte Carlo standard deviation
    intcox_betas_MC_sd <- MC.sd(intcox_betas)
    write.results(results=intcox_betas_MC_sd, data_dir=data_dir, file.name='intcox_betas_MC_sd')
    
    # Compute average asymptotic standard deviation
    intcox_betas_mean_sd <- colMeans(intcox_betas_sd)
    write.results(results=intcox_betas_mean_sd, data_dir=data_dir, file.name='intcox_betas_mean_sd')
    
    # Compute mean squared error
    intcox_betas_mse <- mse(intcox_betas, c(beta1, beta2, beta3))
    write.results(results=intcox_betas_mse, data_dir=data_dir, file.name='intcox_betas_mse')
    
    # Collate MIICD results
    
    # data.frames to accumulate the estimates  
    miicd_betas <- data.frame()
    miicd_betas_sd <- data.frame()
    miicd_surv <- data.frame()
    
    for (i in 1:N) {
      cat(ifelse(i %% 50 == 0, '.\n', '.'))
      
      # process ith miicd replication
      miicd <- get.nth.miicd.result(i, data_dir)
      miicd_betas <- rbind(miicd_betas, miicd$b[,1])
      miicd_betas_sd <- rbind(miicd_betas_sd, miicd$b[,3])
      miicd_surv <- rbind(miicd_surv, miicd$surv[,2])
    }
    colnames(miicd_betas) <- c('beta_1','beta_2','beta_3')
    save(miicd_betas=miicd_betas, file=file.path(data_dir, 'miicd_betas.rda'))
    colnames(miicd_betas_sd) <- c('beta_1','beta_2','beta_3')
    save(miicd_betas_sd=miicd_betas_sd, file=file.path(data_dir, 'miicd_betas_sd.rda'))
    # save survival curve info
    miicd_surv_xs <- miicd$surv[,1]
    colnames(miicd_surv) <- seq_len(ncol(miicd_surv))
    miicd_surv_join <- list(time=miicd_surv_xs, surv=miicd_surv)
    save(miicd_surv_join, file=file.path(data_dir, 'miicd_surv.rda'))
    
    # Compute bias
    miicd_betas_bias <- bias(miicd_betas, c(beta1, beta2, beta3))
    write.results(results=miicd_betas_bias, data_dir=data_dir, file.name='miicd_betas_bias')
    
    # Compute Monte Carlo standard deviation
    miicd_betas_MC_sd <- MC.sd(miicd_betas)
    write.results(results=miicd_betas_MC_sd, data_dir=data_dir, file.name='miicd_betas_MC_sd')
    
    # Compute average asymptotic standard deviation
    miicd_betas_mean_sd <- colMeans(miicd_betas_sd)
    write.results(results=miicd_betas_mean_sd, data_dir=data_dir, file.name='miicd_betas_mean_sd')
    
    # Compute mean squared error
    miicd_betas_mse <- mse(miicd_betas, c(beta1, beta2, beta3))
    write.results(results=miicd_betas_mse, data_dir=data_dir, file.name='miicd_betas_mse')
  }
  
  cat('\nFinished writing files.\n')
}