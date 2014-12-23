# This file constructs a plot of the measured absolute bias in MPL
#  regression coefficient estimates at the three censoring proportions
#  across the three sample sizes

# As the MPL code is not included (currently unpublished), this is intended
# for viewing only

library(ggplot2)
library(reshape)
library(nortest)
options(digits=2, width=128, scipen=15)

# these correspond to the censoring rates (%) and sample sizes
censoring <- c(80, 50, 20)
samplesize <- c(100, 500, 5000)

results <- data.frame()
for (i in 1:9) {
  # the following relies on the naming convention of the simulation result files
  estimate.file <- paste0("run1_00",i,"_results/run1_00",i,'/phmpl_betas.rda')
  load(estimate.file) # loads the data frame phmpl_betas
  
  row = data.frame(pc.censored  = censoring[1 + floor((i-1) / 3)],
                   sample.size  = samplesize[1 + (i-1) %% 3],
                   beta1   = mean(phmpl_betas$beta_1) - 1.5,
                   beta2   = mean(phmpl_betas$beta_2) - 0.5,
                   beta3   = mean(phmpl_betas$beta_3) + 0.1
  )
  results <- rbind(results, row)
}

# bias test results
print(results[,c(1,2,grep("bias.p",colnames(results)))])

# plot
plot.data <- melt(results, id=c(1,2))
plot.data$pc.censored <- factor(with(plot.data, paste0(pc.censored, "% censored")))
ggplot(plot.data, aes(x=sample.size, y=abs(value), color=variable, shape=variable)) +
      geom_point() +
      geom_line() +
      facet_wrap(~pc.censored) +
      scale_x_log10(breaks=samplesize) +
      ggtitle('MPL: absolute bias of regression estimates') +
      ylab('Absolute bias') + xlab('Sample size (log scale)') +
      theme_bw(20) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('asymptotic_bias.pdf', width = 10, height=6)
