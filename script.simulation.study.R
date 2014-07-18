rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/functions/dinvgamma.R')
source('~/functions/rMVN.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')

make.sim.plot <- function(out){
  layout(matrix(1:2, nrow = 1))
  
  matplot(out$fort.raster, type = 'l', ylim = c(min(out$fort.raster) - 2 * max(sqrt(out$var.save)), max(out$fort.raster) + 2 * max(sqrt(out$var.save))), main = "Posterior Predictions", ylab = 'T', xlab = 'location')
  for(i in 1:dim(out$var.save)[2]){
    matplot(out$fort.raster - 2 * sqrt(out$var.save[, i]), type = 'l', add = TRUE, col = adjustcolor(1:10, alpha = 0.25), lty = 'dashed')
    matplot(out$fort.raster + 2 * sqrt(out$var.save[, i]), type = 'l', add = TRUE, col = adjustcolor(1:10, alpha = 0.25), lty = 'dashed')
  }
  
  MSPE <- (out$fort.raster - matrix(unlist(Z.list.hist), nrow = m, byrow = FALSE))^2
  matplot(MSPE, type = 'l', main = paste('MSPE', round(mean(MSPE), 4)), ylab = "MSPE", xlab = "location")
}

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

##
## Plot simulated data
##

jpeg(file = paste('~/1dSpatialSim/plots/SimulationStudy', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
layout(matrix(1:2, nrow = 1))
plot.Y.field(field$Y.list[1:(reps / 2)], field$H.list[1:(reps / 2)], locs, ylab = 'T', xlab = 'location')
plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = 'Full Data', ylab = 'T', xlab = 'location')
dev.off()

Y.list <- field$Y.list[1:(reps / 2)] 
H.list <- field$H.list[1:(reps / 2)] 
Z.list.hist <- field$Z.list[1:(reps / 2)]
Z.list.pca <- field$Z.list[(reps / 2 + 1):reps]
X <- matrix(unlist(Z.list.pca), ncol = reps / 2, byrow = FALSE)
# matplot(X, type = 'l')

num.pca <- 3 


## Simple PCA regression
source('~/1dSpatialSim/pcaRegression/mcmc.pca.R')
##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, num.pca)
sigma.squared.0 <- 25
Sigma.0 <- sigma.squared.0 * diag(num.pca)
alpha.beta <- 100
beta.beta <- 0.1
##
alpha.epsilon <- 0.1
beta.epsilon <- 0.1 
##
alpha.phi <- 0.1 
beta.phi <- 0.1 
##

n.mcmc <- 25000
n.burn <- floor(n.mcmc / 5)

##
## Fit spatial MCMC kriging model
##

out <- mcmc.1d(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)


##
## Plot output
##
jpeg(file = paste('~/1dSpatialSim/plots/pcaRegression', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')

make.sim.plot(out)

dev.off()


##
## Spatial smoothing model
##
source('~/1dSpatialSim/spatialSmoothing/mcmc.spatial.R')

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, dim(X)[2])
Sigma.0 <- sigma.squared.0 * diag(dim(X)[2])
##
alpha.eta <- 0.1
beta.eta <- 0.1
##
sigma.squared.beta.tune <- 0.025
sigma.squared.eta.tune <- 0.25
sigma.squared.epsilon.tune <- 0.5
phi.tune <- 0.25

##
## Fit spatial MCMC kriging model
##

out <- mcmc.1d(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)


jpeg(file = paste('~/1dSpatialSim/plots/spatialRegression', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
make.sim.plot(out)
dev.off()


## 
## Predictive Process
##

source('~/fortTemp/pcaSpatialPPModified/1dmcmc.pca.spatial.pp.R')
##

mu.0 <- rep(0, num.pca)
Sigma.0 <- sigma.squared.0 * diag(num.pca)
D <- as.matrix(dist(locs))
##
## tuning for 99 knot locations
##
s.star <- seq(0.01, 0.99, 0.01)
sigma.squared.eta.tune <- 0.7
sigma.squared.epsilon.tune <- .062575
phi.tune <- 0.375

out <- mcmc.pp(Y.list, H.list, X, locs, D, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, s.star)

jpeg(file = paste('~/1dSpatialSim/plots/PPSimulation99knots', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
make.sim.plot(out)
dev.off()


##
## tuning for 39 knot locations
##
s.star <- seq(0.025, 0.975, 0.025)
sigma.squared.eta.tune <- 0.70
sigma.squared.epsilon.tune <- .070
phi.tune <- 0.6

out <- mcmc.pp(Y.list, H.list, X, locs, D, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, s.star)

jpeg(file = paste('~/1dSpatialSim/plots/PPSimulation39knots', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
make.sim.plot(out)
dev.off()


##
## tuning for 9 knot locations
##
s.star <- seq(0.1, 0.9, 0.1)
sigma.squared.eta.tune <- 0.85
sigma.squared.epsilon.tune <- .085
phi.tune <- 0.35

out <- mcmc.pp(Y.list, H.list, X, locs, D, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, s.star)

jpeg(file = paste('~/1dSpatialSim/plots/PPSimulation9knots', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
make.sim.plot(out)
dev.off()

