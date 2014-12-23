# wblcen: Generate survival data with random independent censoring
#         (weibull distribution)
#
# Inputs: n -- sample size
#         betas -- parameters on covariates
#         k -- shape parameter
#         b -- fixed distance between pairs of left and right interval endpoints
#         epct -- proportion of exact observations
#
# Usage example: wblcen(500, c(1.5,0.5,-0.1), 3, 2, 0.8)

wblcen <- function (n, betas, k, b, epct) {
x1 = runif(n)
x2 = runif(n, max=5)
x3 = runif(n, max=7)
X = cbind(x1, x2, x3)

lambda = (exp(-X %*% matrix(betas)))^(1/k)
t = rweibull(n, shape=k, scale=lambda)
d = rbinom(n, size=1, prob=epct)
cl = runif(n)
cr = cl + b*runif(n)
tl = rep(0, n)
tr = tl

for (i in 1:n){
  if (d[i]==1){
    tl[i] = t[i]
    tr[i] = t[i]
  }
  else if (d[i]==0){
    if (t[i] < cl[i]){
      tl[i] = 0
      tr[i] = cl[i]
    }
    else if (t[i] >= cl[i] && t[i] < cr[i]){
      tl[i] = cl[i]
      tr[i] = cr[i]
    }
    else{
      tl[i] = cr[i]
      tr[i] = Inf
    }
  }
}
sdata = cbind(tl, tr)

return(list(y=sdata, X=X))
}