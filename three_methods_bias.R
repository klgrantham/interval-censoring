# This file constructs a plot of the measured absolute bias in regression
#  coefficient estimates for the MPL, MI and ICM methods.

# As the MPL code is not included (currently unpublished), this is intended
# for viewing only

rm(list=ls())
library(ggplot2)
library(reshape)

samplesize <- c(100, 500, 5000)

# build a row of the result table
res.row <- function(method, n, betas) {
  data.frame(method = method,
             sample.size  = n,
             beta1        = mean(betas$beta_1) - 1.5,
             beta2        = mean(betas$beta_2) - 0.5,
             beta3        = mean(betas$beta_3) + 0.1)
}

results <- data.frame()
for (i in 1:3) {
  # the following relies on the naming convention of the simulation result files
  res.dir <- paste0("complete",i,"_all_results/complete",i,"/")
  load(paste0(res.dir, "phmpl_betas.rda"))
  load(paste0(res.dir, "miicd_betas.rda"))
  load(paste0(res.dir, "intcox_betas.rda"))
  
  n = samplesize[i]
  results <- rbind(results, 
                   res.row("MPL", n, phmpl_betas),
                   res.row("MI", n, miicd_betas),
                   res.row("ICM",n, intcox_betas) )
}

# plot
plot.data <- melt(results, id=c(1,2))
plot.data$method <- with(plot.data, paste(method, 'method'))
ggplot(plot.data, aes(x=sample.size, y=abs(value), colour=variable, shape=variable)) +
  geom_point(size=3) +
  geom_line() +
  facet_wrap(~method) +
  scale_x_log10(breaks=samplesize) +
  ggtitle('Absolute bias of regression estimates') +
  ylab('Absolute bias') + xlab('Sample size (log scale)') +
  theme_bw(20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('three_methods_bias.pdf', width = 10, height=6)
