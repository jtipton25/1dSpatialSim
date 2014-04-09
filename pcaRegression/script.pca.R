rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.pca.R')
# source('make.spatial.field.pca.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/pcaRegression/mcmc.pca.R')

##
## Simulate Data
##

##
## Plot of 1d spatial mcmc output
##

make.output.plot <- function(out){
  n.burn <- floor(n.mcmc / 10)
  #x11()
  layout(matrix(1:16, nrow = 4))
  #
  matplot(t(out$mu.beta.save)[(n.burn + 1):n.mcmc, ], type = 'l')
  abline(h = beta[1], col = 'black')
  abline(h = beta[2], col = 'red')
  #
  plot(out$sigma.squared.beta.save[(n.burn + 1):n.mcmc], type = 'l')
  #
  plot(out$sigma.squared.epsilon.save[(n.burn + 1):n.mcmc], type = 'l')
  abline(h = s2.e, col = 'red')
  #
  matplot(out$fort.raster, type = 'l')
  #
  plot.Z.field(Z.list.hist, locs = locs, main = "True Surface")
  #
  plot.Z.field(Z.list.pca, locs = locs, main = "Pallette of Signals")
  #
  plot.Y.field(Y.list, H.list, locs = locs)
  #
  hist(out$mu.beta.save[1, ][(n.burn + 1):n.mcmc])
  abline(v = beta[1], col = 'red')
  abline(v = quantile(out$mu.beta.save[1, ], probs = c(0.025, 0.975)), col = 'blue')
  #
  hist(out$mu.beta.save[2, ][(n.burn + 1):n.mcmc])
  abline(v = beta[2], col = 'red')
  abline(v = quantile(out$mu.beta.save[2, ], probs = c(0.025, 0.975)), col = 'blue')
  #
  MSPE <- (out$fort.raster - matrix(unlist(Z.list.hist), nrow = m, byrow = FALSE))^2
  matplot(MSPE, type = 'l', main = 'MSPE')
}




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
plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = 'Full Data')

Y.list <- field$Y.list[1:(reps / 2)] 
H.list <- field$H.list[1:(reps / 2)] 
Z.list.hist <- field$Z.list[1:(reps / 2)]
Z.list.pca <- field$Z.list[(reps / 2 + 1):reps]
X <- matrix(unlist(Z.list.pca), ncol = reps / 2, byrow = FALSE)
matplot(X, type = 'l')
num.pca <- 3 

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, num.pca)#rep(0, dim(X)[2])
sigma.squared.0 <- 25
Sigma.0 <- sigma.squared.0 * diag(num.pca)
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
sigma.squared.eta.tune <- 0.5
sigma.squared.epsilon.tune <- 0.275
phi.tune <- 0.5

n.mcmc <- 5000

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
# Rprof(file = 'spatial.Rprof.out')
out.pca <- mcmc.1d(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
# Rprof(NULL)
# summaryRprof('spatial.Rprof.out')
finish <- Sys.time() - start
finish 

#5000 iterations takes 12.15 minutes for m = 1000 and reps = 10
 
##
## Plot output
##

make.output.plot(out.pca)

apply(out.pca$mu.beta.save[, (n.mcmc / 10 + 1):n.mcmc], 1, mean)

make.MSPE <- function(s, out){
  (out$fort.raster[, s] - Z.list.hist[[s]])^2
}
MSPE.pca <- sapply(1:(reps /  2), make.MSPE, out = out.pca)
# matplot(MSPE.pca, type = 'l', main = 'MSPE')
mean(MSPE.pca)
