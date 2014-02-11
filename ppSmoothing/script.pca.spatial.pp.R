rm(list = ls())
set.seed(1)

##
## Libraries and Subroutines
##
setwd('~/1dSpatialSim/')
source('dinvgamma.R')
library(mvtnorm)
source('make.output.plot.R')
source('make.spatial.field.pca.R')
setwd('~/1dSpatialSim/ppSmoothing//')
source('mcmc.pca.spatial.pp.R')
# source('mcmc.spatial.pp.R')

##
## Simulate Data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 200 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.01
samp.size <- 5:40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)
plot.Y.field(field$Y.list, field$H.list, locs)
plot.Z.field(field$Z.list, locs)

Y.list <- field$Y.list[1:100]
H.list <- field$H.list[1:100]
Z.list <- field$Z.list[101:200]

matplot(matrix(unlist(Z.list), ncol = 100, byrow = FALSE), type = 'l')
X <- matrix(unlist(Z.list), ncol = 100, byrow = FALSE)
num.pca <- 3

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, 3)#c(0, 2)#rep(0, dim(X)[2])
sigma.squared.0 <- 0.025
Sigma.0 <- sigma.squared.0 * diag(num.pca)
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
sigma.squared.eta.tune <- 0.0075
sigma.squared.epsilon.tune <- 0.005
phi.tune <- 0.75

##
## Knots for predictive process
##

s.star <- seq(0.1, 0.9, 0.1)

sigma.squared.eta.tune <- 0.00275
sigma.squared.epsilon.tune <- 0.0020
phi.tune <- 0.0050

n.mcmc <- 1000

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

make.output.plot(out)
apply(out$mu.beta.save, 1, mean)
