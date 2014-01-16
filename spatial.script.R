rm(list = ls())
set.seed(203)

##
## Libraries and Subroutines
##

dinvgamma = function(x, shape = 1, rate = 1, scale = 1/rate, log = FALSE) {
  # return( rate^shape / gamma(shape) * exp( - rate / x) * x^( - shape - 1))
  logval = shape * log(rate) - lgamma(shape) - rate / x - (shape + 1) * log(x)
  if (log)
    return(logval)
  else
    return(exp(logval))
}

##
## Simulate Data
##

m <- 100 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 100 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 2
s2.e <- 0.5
samp.size <- 40

source('make.spatial.field.R')
field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)


##
## Initialize priors and tuning paramteters
##

mu.0 <- c(0, 2)#rep(0, dim(X)[2])
sigma.squared.0 <- 0.025
#Sigma.0 <-
alpha.beta <- 10
beta.beta <- 10
curve(dinvgamma(x, alpha.beta, beta.beta))
alpha.eta <- 10
beta.eta <- 10
curve(dinvgamma(x, alpha.eta, beta.eta))
alpha.epsilon <- 40
beta.epsilon <- 20
curve(dinvgamma(x, alpha.epsilon, beta.epsilon))
alpha.phi <- 10
beta.phi <- 10
curve(dinvgamma(x, alpha.phi, beta.phi))

sigma.squared.beta.tune <- 0.025
sigma.squared.eta.tune <- 0.125
sigma.squared.epsilon.tune <- 0.075
phi.tune <- 0.25

n.mcmc <- 1000

source('mcmc.spatial.R')

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
finish <- Sys.time() - start
finish #100 iterations takes 8 minutes

##
## Plot output
##

#x11()
layout(matrix(1:9, nrow = 3))
matplot(t(out$mu.beta.save), type = 'l')
abline(h = beta[1], col = 'black')
abline(h = beta[2], col = 'red')
plot(out$sigma.squared.beta.save, type = 'l')
plot(out$sigma.squared.epsilon.save, type = 'l', main = paste("accept rate", round(out$epsilon.accept, 2)))
abline(h = s2.e)
plot(out$sigma.squared.eta.save, type = 'l', main = paste("accept rate", round(out$eta.accept, 2)))
abline(h = s2.s)
plot(out$phi.save, type = 'l', main = paste("accept rate", round(out$phi.accept, 2)))
abline(h = phi)
matplot(out$fort.raster, type = 'l')
plot.field(field$Z.list, H.list = rep(list(1:length(field$Z.list[[1]])), reps), locs = locs)
#plot.field(Z.list, H.list = rep(list(1:length(Z.list[[1]])), reps), locs = locs)
hist(out$mu.beta.save[1, ])
abline(v = mean(out$mu.beta.save[1, ]), col = 'red')
abline(v = quantile(out$mu.beta.save[1, ], probs = c(0.025, 0.975)), col = 'blue')
hist(out$mu.beta.save[2, ])
abline(v = mean(out$mu.beta.save[2, ]), col = 'red')
abline(v = quantile(out$mu.beta.save[2, ], probs = c(0.025, 0.975)), col = 'blue')

