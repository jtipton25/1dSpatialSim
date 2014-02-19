rm(list = ls())
set.seed(10)

##
## Libraries and Subroutines
##
setwd('~/1dSpatialSim/')
source('dinvgamma.R')
# source('make.output.plot.R')
source('make.output.plot.ci.R')
library(mvtnorm)
source('make.spatial.field.R')
setwd('~/1dSpatialSim/spatialSmoothingFixed//')
source('mcmc.spatial.R')


##
## Plot of 1d spatial mcmc output
##

make.output.plot <- function(out){
  n.burn <- floor(n.mcmc / 5)
  #x11()
  layout(matrix(1:16, nrow = 4))
  #
#   matplot(t(out$mu.beta.save)[(n.burn + 1):n.mcmc, ], type = 'l')
#   abline(h = beta[1], col = 'black')
#   abline(h = beta[2], col = 'red')
  #
  plot(out$sigma.squared.beta.save[(n.burn + 1):n.mcmc], type = 'l')
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
  matplot(out$fort.raster, type = 'l')#, ylim = c(min(out$fort.raster) - 2*max(sqrt(out$var.save)), max(out$fort.raster) + 2*max(sqrt(out$var.save))))
#   matplot(out$fort.raster - 2*sqrt(out$var.save), type = 'l', add = TRUE, col = 'red', lty = 'dashed')
#   matplot(out$fort.raster + 2*sqrt(out$var.save), type = 'l', add = TRUE, col = 'red', lty = 'dashed')
  points(X %*% beta, type = 'l', col = 'red')  
  #
  plot.Z.field(field$Z.list, locs = locs, main = "True Surface")
  #
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

m <- 500 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 20 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.01
samp.size <- 5:40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

layout(matrix(1:2, ncol = 2))
# plot.Y.field(field$Y.list[1:(reps / 2)], field$H.list[1:(reps / 2)], locs)
# plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = "Full Data")
plot.Y.field(field$Y.list, field$H.list, locs)
plot.Z.field(field$Z.list, locs, main = "Full Data")

# Y.list <- field$Y.list[1:(reps / 2)]
# H.list <- field$H.list[1:(reps / 2)]
# Z.list <- field$Z.list[(reps / 2 + 1):reps]
Y.list <- field$Y.list
H.list <- field$H.list
Z.list <- field$Z.list

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, dim(X)[2])
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

n.mcmc <- 5000

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
 x11()
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
