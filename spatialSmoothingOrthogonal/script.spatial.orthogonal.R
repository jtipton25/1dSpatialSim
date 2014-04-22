rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
setwd('~/1dSpatialSim/spatialSmoothingOrthogonal/')
source('~/1dSpatialSim/spatialSmoothingOrthogonal/mcmc.spatial.orthogonalEdit.R')

##
## Simulate Data
##

m <- 250 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 20 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.01
samp.size <- 5:40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

Y.list <- field$Y.list[1:(reps / 2)] 
H.list <- field$H.list[1:(reps / 2)] 
Z.list.hist <- field$Z.list[1:(reps / 2)]
Z.list.pca <- field$Z.list[(reps / 2 + 1):reps]
X <- matrix(unlist(Z.list.pca), ncol = reps / 2, byrow = FALSE)
# matplot(X, type = 'l')
num.pca <- 3 

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, num.pca)
sigma.squared.0 <- 25
Sigma.0 <- sigma.squared.0 * diag(num.pca)
alpha.beta <- 20
beta.beta <- 0.2
#curve(dinvgamma(x, alpha.beta, beta.beta))
##
alpha.eta <- 12
beta.eta <- 12
#curve(dinvgamma(x, alpha.eta, beta.eta), from = 0, to = 6)
#abline(v = s2.s, col = 'red')
##
alpha.epsilon <- 3
beta.epsilon <- 2
#curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 6)
#abline(v = s2.e, col = 'red')
##
alpha.phi <- 10
beta.phi <- 20
#curve(dinvgamma(x, alpha.phi, beta.phi), from = 0, to = 6)
#abline(v = phi, col = 'red')
##
sigma.squared.eta.tune <- 0.000075
sigma.squared.epsilon.tune <- 0.0005
phi.tune <- 0.0050

n.mcmc <- 5000

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
# Rprof(filename = '~/1dSpatialSim/spatialOrthogonalRprof.out', line.profiling = TRUE)
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
# Rprof(NULL)
# summaryRprof(filename = '~/1dSpatialSim/spatialOrthogonalRprof.out', lines = 'show')

finish <- Sys.time() - start
finish #500 iterations takes 2.23 minutes for m = 100 and reps = 100
#500 iterations takes 5.3 minutes for m = 1000 and reps = 100

##
## Plot output
##

make.output.plot(out)

## identifiability between beta_0 and sigma.squared.epsilon???
#matplot(out$beta.save[1, , (n.mcmc / 10 + 1):n.mcmc], type = 'l', ylim = c(min(out$beta.save[, , (n.mcmc / 10 + 1):n.mcmc]), max(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc])))
#matplot(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc], type = 'l', add = TRUE)

#apply(out$mu.beta.save[, (n.mcmc / 10 + 1):n.mcmc], 1, mean)
#matplot(out$fort.raster, type = 'l', main = 'Posterior Predictive')

#MSPE <- (out$fort.raster - matrix(unlist(field$Z.list), nrow = m, byrow = FALSE))^2
#matplot(MSPE, type = 'l', main = 'MSPE')

