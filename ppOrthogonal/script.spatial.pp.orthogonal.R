rm(list = ls())
set.seed(203)

##
## Libraries and Subroutines
##

setwd('~/1dSpatialSim/')
source('dinvgamma.R')
source('rMVN.R')
#library(mvtnorm)
source('make.output.plot.R')
source('make.spatial.field.R')
setwd('~/1dSpatialSim/ppOrthogonal/')
source('mcmc.spatial.pp.orthogonal.R')

##
## Simulate Data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 100 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 2
s2.e <- 0.5
samp.size <- 40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

##
## Initialize priors and tuning paramteters
##

mu.0 <- c(0, 2)#rep(0, dim(X)[2])
sigma.squared.0 <- 0.025
#Sigma.0 <- sigma.squared.0 * diag(dim(X)[2])
alpha.beta <- 2
beta.beta <- 0.2
curve(dinvgamma(x, alpha.beta, beta.beta))
##
alpha.eta <- 12
beta.eta <- 12
curve(dinvgamma(x, alpha.eta, beta.eta), from = 0, to = 6)
abline(v = s2.s, col = 'red')
##
alpha.epsilon <- 3
beta.epsilon <- 2
curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 6)
abline(v = s2.e, col = 'red')
##
alpha.phi <- 10
beta.phi <- 20
curve(dinvgamma(x, alpha.phi, beta.phi), from = 0, to = 6)
abline(v = phi, col = 'red')
##
sigma.squared.eta.tune <- 0.00075
sigma.squared.epsilon.tune <- 0.000025
phi.tune <- 0.00075

n.mcmc <- 1000

##
## Knots for predictive process
##

s.star <- seq(0.1, 0.9, 0.1)

##
## Fit spatial MCMC kriging model
##

## Something is happening to sigma.squared.epsilon.save that is driving it to 0...

##
## Profile the code
##
Rprof(filename = 'Rprof.pp.out', line.profiling = TRUE)
start <- Sys.time()
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, s.star)
finish <- Sys.time() - start
finish 
Rprof(NULL)

summaryRprof(filename = 'Rprof.pp.out', lines = 'show')

#500 iterations takes 2.21 minutes for m = 100 and reps = 100
#500 iterations takes 3.6 minutes for m = 1000 and reps = 100

##
## Plot output
##

make.output.plot(out)
matplot(out$beta.save[1, , (n.mcmc / 10 + 1):n.mcmc], type = 'l', ylim = c(min(out$beta.save[, , (n.mcmc / 10 + 1):n.mcmc]), max(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc])))
matplot(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc], type = 'l', add = TRUE)

apply(out$mu.beta.save[, (n.mcmc / 10 + 1):n.mcmc], 1, mean)
matplot(out$fort.raster, type = 'l', main = 'Posterior Predictive')

MSPE <- (out$fort.raster - matrix(unlist(field$Z.list), nrow = m, byrow = FALSE))^2
matplot(MSPE, type = 'l', main = 'MSPE')
out$accept
