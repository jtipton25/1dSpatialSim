rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/spatialSmoothingPP/mcmc.spatial.pp.R')

##
## Simulate Data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 20 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.1
samp.size <- 5:40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

layout(matrix(1:2, ncol = 2))
plot.Y.field(field$Y.list[1:(reps / 2)], field$H.list[1:(reps / 2)], locs)
plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = "Full Data")

Y.list <- field$Y.list[1:(reps / 2)]
H.list <- field$H.list[1:(reps / 2)]
Z.list.hist <- field$Z.list[1:(reps / 2)]
Z.list.pca <- field$Z.list[(reps / 2 + 1):reps]
# X <- matrix(unlist(Z.list), ncol = reps / 2, byrow = FALSE)
# matplot(X, type = 'l')

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, dim(X)[2])
sigma.squared.0 <- 25
Sigma.0 <- sigma.squared.0 * diag(dim(X)[2])
alpha.beta <- 2
beta.beta <- 0.2
curve(dinvgamma(x, alpha.beta, beta.beta))
##
alpha.eta <- 12
beta.eta <- 12
curve(dinvgamma(x, alpha.eta, beta.eta), from = 0, to = 6)
abline(v = s2.s, col = 'red')
##
alpha.epsilon <- 1
beta.epsilon <- 0.1
curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 6)
abline(v = s2.e, col = 'red')
##
alpha.phi <- 4
beta.phi <- 1
curve(dinvgamma(x, alpha.phi, beta.phi), from = 0, to = 6)
abline(v = phi, col = 'red')
##

##
## Knots for predictive process
##

s.star <- seq(0.1, 0.9, 0.1)

sigma.squared.eta.tune <- 0.275
sigma.squared.epsilon.tune <- 0.0250
phi.tune <- 75.50

n.mcmc <- 20000

##
## Fit spatial MCMC kriging model
##

## Not sure I believe the posterior predictive values...

start <- Sys.time()
# Rprof(file = 'spatial.pp.Rprof.out')
out <- mcmc.1d(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, s.star)
# Rprof(NULL)
# summaryRprof('spatial.pp.Rprof.out')
finish <- Sys.time() - start
finish 

## 5000 iterations takes 4.45 minutes for predictive process for m = 1000 and reps = 10

##
## Plot output
##
# x11()
make.output.plot(out)


apply(out$mu.beta.save, 1, mean)
str(out$mu.beta.save)

MSPE <- matrix(nrow = length(locs), ncol = reps / 2)
for(i in 1:(reps / 2)){
  MSPE[, i] <- (out$fort.raster[, i] - Z.list.hist[[i]])^2
}
matplot(MSPE, type = 'l')
mean(MSPE)

