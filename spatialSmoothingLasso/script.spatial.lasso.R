rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##
source('~/1dSpatialSim/functions/dinvgamma.R')
# source('make.output.plot.R')
library(mvtnorm)
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/spatialSmoothingLassoFixedEffect/mcmc.spatial.lasso.R')

##
## Plot of 1d spatial mcmc output
##

make.output.plot <- function(out){
  n.burn <- floor(n.mcmc / 10)
  layout(matrix(1:9, nrow = 3))
  #
  plot(out$sigma.squared.epsilon.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$epsilon.accept, 2)))
  abline(h = s2.e, col = 'red')
  #
  plot(out$sigma.squared.eta.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$eta.accept, 2)))
  abline(h = s2.s, col = 'red')
  #
  plot(out$phi.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$phi.accept, 2)))
  abline(h = phi, col = 'red')
  #
  matplot(out$fort.raster, type = 'l')
  #
  plot.Z.field(field$Z.list, locs = locs, main = "True Surface")
  #
  plot.Y.field(field$Y.list, field$H.list, locs = locs)
  #
  hist(out$beta.save[1, , ][(n.burn + 1):n.mcmc])
  abline(v = beta[1], col = 'red')
  abline(v = quantile(out$beta.save[1, , ], probs = c(0.025, 0.975)), col = 'blue')
  #
  hist(out$beta.save[2, , ][(n.burn + 1):n.mcmc])
  abline(v = beta[2], col = 'red')
  abline(v = quantile(out$beta.save[2, , ], probs = c(0.025, 0.975)), col = 'blue')
  #
  MSPE <- (out$fort.raster - matrix(unlist(field$Z.list), nrow = m, byrow = FALSE))^2
  matplot(MSPE, type = 'l', main = 'MSPE')
}


##
## Simulate Data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 2 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.1
samp.size <- 40:100

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

layout(matrix(1:2, ncol = 2))
plot.Z.field(field$Z.list, locs, main = "Actual data")
plot.Y.field(field$Y.list, field$H.list, locs)
Y.list <- field$Y.list
H.list <- field$H.list
##
## Initialize priors and tuning paramteters
##

mu.beta <- beta
mu.0 <- beta #rep(0, dim(X)[2])
sigma.squared.0 <- 0.025
#Sigma.0 <-
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


alpha.lambda <- 30#1 
beta.lambda <- 1#30
curve(dgamma(x, alpha.lambda, beta.lambda))
# alpha.lambda = 1, beta.lambda = 30 provides little shrinkage
# alpha.lambda = 30, beta.lambda = 1 provides strong shrinkage

##
sigma.squared.beta.tune <- 0.05
sigma.squared.eta.tune <- 0.5
sigma.squared.epsilon.tune <- 0.0150
phi.tune <- 0.75
sigma.squared.gamma.tune <- 0.3
sigma.squared.gamma.tune <- 0.0015
n.mcmc <- 5000

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, n.mcmc, mu.beta, alpha.epsilon, beta.epsilon, alpha.phi, beta.phi, alpha.lambda, beta.lambda, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, sigma.squared.gamma.tune)
finish <- Sys.time() - start
finish 

#5000 iterations takes 12.15 minutes for m = 1000 and reps = 10

##
## Plot output
##
x11()
make.output.plot(out)
out$gamma.accept
## identifiability between beta_0 and sigma.squared.epsilon???
# matplot(out$beta.save[1, , (n.mcmc / 10 + 1):n.mcmc], type = 'l', ylim = c(min(out$beta.save[, , (n.mcmc / 10 + 1):n.mcmc]), max(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc])))
# matplot(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc], type = 'l', add = TRUE)

# apply(out$mu.beta.save[, (n.mcmc / 10 + 1):n.mcmc], 1, mean)
# matplot(out$fort.raster, type = 'l', main = 'Posterior Predictive')

# MSPE <- (out$fort.raster - matrix(unlist(field$Z.list), nrow = m, byrow = FALSE))^2
# matplot(MSPE, type = 'l', main = 'MSPE')
