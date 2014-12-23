# This file contains the following functions:
#   plot.preprocess -- preprocesses the output files from any method
#     and creates a data frame that can be used for plotting
#   plot.generate.MPL.basehaz -- generates a baseline hazard plot
#     for the MPL method (MPL is the only method that estimates a smooth baseline
#     hazard) with both asymptotic and Monte Carlo simultaneous confidence intervals
#   plot.generate.MI.ICM.cumbasehaz -- generates a cumulative baseline hazard
#     plot for the MI and ICM methods with Monte Carlo simultaneous confidence intervals

# Note: Code for the MPL method is not yet publicly available, so plots for
#        this method cannot be generated at this time.
generate.MPL.plots <- FALSE

# function to preprocess output files and return data frame for plotting
plot.preprocess <- function (time, basehaz, basehaz_sd=NULL) {
  N <- dim(basehaz)[1] # Number of replications for phmpl
  M <- dim(basehaz)[2] # Number of points at which curve has been computed
  # Monte Carlo simultaneous confidence intervals
  basehaz_mins <- matrix(0, 1, M)
  basehaz_maxs <- matrix(0, 1, M)
  avg_basehaz <- matrix(0, 1, M)
  for (j in 1:M) {
    basehaz_locj <- sort(basehaz[,j])
    basehaz_mins[,j] <- basehaz_locj[ceiling(N*0.025)]
    basehaz_maxs[,j] <- basehaz_locj[floor(N*0.975)]
    basehaz_mid95 <- basehaz_locj[ceiling(N*0.025):floor(N*0.975)]
    avg_basehaz[,j] <- mean(basehaz_mid95, na.rm=TRUE)
  }
  # Asymptotic simultaneous confidence intervals
  lbs_asymp <- rep(NA, N)
  ubs_asymp <- rep(NA, N)
  if (!is.null(basehaz_sd)) {
    lbs <- matrix(0, dim(basehaz)[1], M)
    ubs <- matrix(0, dim(basehaz)[1], M)
    for (j in 1:M) {
      lbs[,j] <- basehaz[,j] - 1.96*basehaz_sd[,j]
      ubs[,j] <- basehaz[,j] + 1.96*basehaz_sd[,j]
    }
    lbs_asymp <- colMeans(lbs, na.rm=TRUE)
    ubs_asymp <- colMeans(ubs, na.rm=TRUE)
  }
  # Create data frame for plotting
  true_basehaz <- 3*time^2
  true_cumbasehaz <- time^3
  df <- data.frame(x=time, y=as.numeric(avg_basehaz),
                   ytruebh=true_basehaz, ytruebch=true_cumbasehaz,
                   lbs_MC=as.numeric(basehaz_mins), ubs_MC=as.numeric(basehaz_maxs),
                   lbs_asymp=as.numeric(lbs_asymp), ubs_asymp=as.numeric(ubs_asymp))
  return(df)
}

# function to generate baseline hazard plot for MPL method
# code is given but plots are not generated due to 
plot.generate.MPL.basehaz <- function (MPL_df) {
  suppressMessages(require(ggplot2))
  p <- ggplot(MPL_df, aes(x)) +
    geom_line(aes(y=ytruebh), color="black") +
    geom_line(aes(y=y), color="red") +
    geom_line(aes(y=lbs_MC), linetype="dashed", color="red") +
    geom_line(aes(y=ubs_MC), linetype="dashed", color="red") +
    xlab("time") +
    ylab("hazard rate") +
    ggtitle("MPL: baseline hazard estimates") +
    theme_bw(20)
  if (with(MPL_df, all(!is.na(c(lbs_asymp, ubs_asymp))))) {
    p <- p + geom_line(aes(y=lbs_asymp), linetype="longdash", color="purple") +
      geom_line(aes(y=ubs_asymp), linetype="longdash", color="purple")
  }
  print(p)
}

# function to generate cumulative baseline hazard plots for MI and ICM methods
plot.generate.MI.ICM.cumbasehaz <- function (MI_df, ICM_df) {
  suppressMessages(require(ggplot2))
  df <- MI_df
  df$y_ICM <- ICM_df$y
  df$lbs_MC_ICM <- ICM_df$lbs_MC
  df$ubs_MC_ICM <- ICM_df$ubs_MC
  p <- ggplot(df, aes(x)) +
    geom_line(aes(y=ytruebch), color="black") +
    geom_line(aes(y=y), color="blue") +
    geom_line(aes(y=lbs_MC), linetype="dashed", color="blue") +
    geom_line(aes(y=ubs_MC), linetype="dashed", color="blue") +
    geom_line(aes(y=y_ICM), color="yellow") +
    geom_line(aes(y=lbs_MC_ICM), linetype="dashed", color="yellow") +
    geom_line(aes(y=ubs_MC_ICM), linetype="dashed", color="yellow") +
    xlab("time") +
    ylab("cumulative hazard rate") +
    ggtitle("Baseline cumulative hazards") +
    theme_bw(20)
  print(p)
}

# load in MPL output files if available and generate plot
if (generate.MPL.plots==TRUE) {
  load('phmpl_bh.rda')
  load('phmpl_bh_sd.rda')
  mpl_bh_time <- phmpl_bh$time
  mpl_bh <- as.matrix(phmpl_bh$bh)
  mpl_bh_sd <- as.matrix(phmpl_bhmod_sd)

  load('phmpl_bch.rda')
  load('phmpl_bch_sd.rda')
  mpl_bch_time <- phmpl_bch$time
  mpl_bch <- as.matrix(phmpl_bch$bch)
  mpl_bch_sd <- as.matrix(phmpl_bchmod_sd)
  
  MPL_df <- plot.preprocess(mpl_bh_time, mpl_bh, mpl_bh_sd)
  plot.generate.MPL.basehaz(MPL_df)
}

# load in MI output files
load('miicd_surv.rda')
miicd_time <- miicd_surv_join$time
miicd_surv <- as.matrix(miicd_surv_join$surv)
miicd_bch <- -log(miicd_surv)

# load in ICM output files
load('intcox_bch.rda')
intcox_time <- intcox_bch_join$time
intcox_bch <- as.matrix(intcox_bch_join$cbh)

# plot cumulative baseline hazards: MI and ICM methods
MI_df <- plot.preprocess(miicd_time, miicd_bch)
ICM_df <- plot.preprocess(intcox_time, intcox_bch)
plot.generate.MI.ICM.cumbasehaz(MI_df, ICM_df)
