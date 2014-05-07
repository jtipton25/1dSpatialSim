rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
# source('make.spatial.field.pca.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/pcaSpatialSmoothing/mcmc.pca.spatial.R')

##
## Simulate Data
##

m <- 20000#1000 # number of spatial locations
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
curve(dinvgamma(x, alpha.beta, beta.beta))
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
sigma.squared.epsilon.tune <- 0.0575
phi.tune <- 0.25

# n.mcmc <- 20000
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

#5000 iterations takes 12.15 minutes for m = 20000 and reps = 20

##
## Plot output
##

# jpeg(file = '~/1dSpatialSim/plots/pcaRegressionSpatial_4_11_2014.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')
make.output.plot(out.pca)
# dev.off()


