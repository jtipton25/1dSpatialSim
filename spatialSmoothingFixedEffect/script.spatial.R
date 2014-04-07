##
## This model is severely overestimating the nugget variance. This is leading to poor predictions...
##

rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/spatialSmoothingFixedEffect/mcmc.spatial.R')

##
## Simulate Data
##

m <- 500 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 2 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.01
samp.size <- 5:40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

layout(matrix(1:2, ncol = 2))
plot.Y.field(field$Y.list, field$H.list, locs)
plot.Z.field(field$Z.list, locs, main = "Full Data")

Y.list <- field$Y.list
H.list <- field$H.list
Z.list <- field$Z.list

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, dim(X)[2])
sigma.squared.0 <- 25
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
sigma.squared.beta.tune <- 0.025
sigma.squared.eta.tune <- 0.25
sigma.squared.epsilon.tune <- 0.075
phi.tune <- 0.25

n.mcmc <- 60

mu.beta <- c(0, 2)

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, mu.beta, n.mcmc, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
finish <- Sys.time() - start
finish 

#5000 iterations takes 12.15 minutes for m = 1000 and reps = 10

##
## Plot output
##
#  x11()
make.output.plot(out)


hist(out$beta.save[1, , ])
abline(v = beta[1], col = 'red')
hist(out$beta.save[2, , ])
abline(v = beta[2], col = 'red')

dim(out$beta.save)
beta.model <- apply(out$beta.save, 1, mean)
phi.model <- mean(out$phi.save)
sigma.squared.epsilon.model <- mean(out$sigma.squared.epsilon.save)
sigma.squared.eta.model <- mean(out$sigma.squared.eta.save)
D <- as.matrix(dist(locs))
model.fit <- t(X %*% as.matrix(beta.model)) + rmvnorm(1, rep(0, m), sigma.squared.eta.model * exp( - D / phi.model))+ rnorm(m, 0, sigma.squared.epsilon.model)
plot(model.fit[1,] ~ locs, type = 'l')
points(locs, X %*% beta, type = 'l', col = 'red')

