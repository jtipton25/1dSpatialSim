rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/pcaSpatialPPModifiedNew/mcmc.pca.spatial.pp.R')

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

# layout(matrix(1:2, ncol = 2))
# plot.Y.field(field$Y.list[1:(reps / 2)], field$H.list[1:(reps / 2)], locs)
# plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = "Full Data")

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

mu.0 <- rep(0, num.pca)#c(0, 2)#rep(0, dim(X)[2])
sigma.squared.0 <- 25
Sigma.0 <- sigma.squared.0 * diag(num.pca)
#Sigma.0 <- sigma.squared.0 * diag(dim(X)[2])
alpha.beta <- 20
beta.beta <- 0.2
##
alpha.eta <- 12
beta.eta <- 12
##
alpha.epsilon <- 3
beta.epsilon <- 2
##
alpha.phi <- 10
beta.phi <- 20
##

##
## Knots for predictive process
##

# s.star <- seq(0.025, 0.975, 0.025)

# sigma.squared.eta.tune <- 0.025
# sigma.squared.epsilon.tune <- 0.00250
# phi.tune <- 0.025

##
## tuning for 99 knot locations
##
s.star <- seq(0.01, 0.99, 0.01)
sigma.squared.eta.tune <- 0.7
sigma.squared.epsilon.tune <- .062575
phi.tune <- 0.375

##
## tuning for 39 knot locations
##
# s.star <- seq(0.025, 0.975, 0.025)
# sigma.squared.eta.tune <- 0.70
# sigma.squared.epsilon.tune <- .070
# phi.tune <- 0.6

##
## tuning for 9 knot locations
##
# s.star <- seq(0.1, 0.9, 0.1)
# sigma.squared.eta.tune <- 0.85
# sigma.squared.epsilon.tune <- .085
# phi.tune <- 0.35

n.mcmc <-20000

##
## Fit spatial MCMC kriging model
##

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

# jpeg(file = '~/1dSpatialSim/plots/pcaRegressionSpatialPPModified99knots_4_14_2014.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')
# jpeg(file = '~/1dSpatialSim/plots/pcaRegressionSpatialPPModified39knots_4_14_2014.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')
# jpeg(file = '~/1dSpatialSim/plots/pcaRegressionSpatialPPModified9knots_4_14_2014.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')

make.output.plot(out)
# dev.off()

# apply(out$mu.beta.save, 1, mean)
# 
# mean(out$phi.save[(n.mcmc / 5 + 1):n.mcmc])
# phi
# mean(out$sigma.squared.epsilon.save[(n.mcmc / 5 + 1):n.mcmc])
# s2.e
# mean(out$sigma.squared.eta.save[(n.mcmc / 5 + 1):n.mcmc])
# s2.s
# MSPE <- (out$fort.raster - matrix(unlist(Z.list.hist), nrow = m))^2
# matplot(MSPE, type = 'l', main = paste(round(mean(MSPE), 4))) ## prediction error
