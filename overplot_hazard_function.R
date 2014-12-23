# This file constructs a plot of estimated baseline hazard functions for
#  the MPL method. It uses a 10% sample of the 1,000 replications.

# As the MPL code is not included (currently unpublished), this is intended
# for viewing only

rm(list=ls())
library(ggplot2)
library(reshape)
load('phmpl_ic_ms_001_cube.rda') # all.dat is a 100x2x1000 array

fraction = 0.10
N=dim(all.dat)[3]
sample.dat <- all.dat[,,sample(1:N,size=fraction*N)]

# compute time range of intersection
min.max.t <- 100
max.min.t <- 0
for (i in 1:N) {
  # left boundary
  max.min.t <- max(max.min.t, Re(all.dat[1,'x',i]))
  # right boundary
  min.max.t <- min(min.max.t, Re(all.dat[100,'x',i]))
}

all.long <- cbind(melt(sample.dat[,1,]), melt(sample.dat[,2,])[,3])
colnames(all.long) <- c('periods','reps','t','bh')

all.long$t <- abs(all.long$t)
all.long$bh <- abs(all.long$bh)
all.long$reps <- factor(all.long$reps)
ggplot(all.long, aes(x=t, y=bh, group=reps, colour=reps)) +
  geom_path(alpha=0.5) +
  geom_vline(x=c(max.min.t,min.max.t),linetype='dashed') +
  ggtitle(expression(atop('MPL: estimated baseline hazard functions',
                          scriptstyle('10% sample of 1,000 replications')))) +
  xlab('time') + ylab('baseline hazard') +
  annotate('text', label='Largest common\nobserved time', x=1.7, y=24) +
  theme_bw(20) +
  theme(legend.position='none')
ggsave('sample_bh_plot.pdf', width=10, height=6)
