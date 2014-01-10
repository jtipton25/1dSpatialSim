
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



params = c(mu.0, Sigma.0, alpha.beta, beta.beta, alpha.eta, beta.eta, alpha.epsilon, beta.epsilon, alpha.phi, beta.phi, n.mcmc = 5000)
mu.0 <- rep(0, dim(X)[2])
sigma.squared.0 <- 100
Sigma.0 <-
alpha.beta <- 10
beta.beta <- 10
curve(dinvgamma(x, alpha.beta, beta.beta))
alpha.eta <- 10
beta.eta <- 10
curve(dinvgamma(x, alpha.eta, beta.eta))
alpha.epsilon <- 10
beta.epsilon <- 10
curve(dinvgamma(x, alpha.epsilon, beta.epsilon))
alpha.phi <- 10
beta.phi <- 10
curve(dinvgamma(x, alpha.phi, beta.phi))
n.mcmc <- 100

sigma.squared.beta.tune <- 0.025
sigma.squared.eta.tune <- 0.025
sigma.squared.epsilon.tune <- 0.025
phi.tune <- 0.25

source('mcmc.spatial.R')

start <-Sys.time()
out <- mcmc.1d(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, sigma.squared.beta, sigma.squared.eta, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta,  s.star, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
finish <- Sys.time() - start
finish #100 iterations takes 8 minutes

layout(matrix(1:8, nrow = 2))
plot(out$sigma.squared.beta.save, type = 'l')
plot(out$sigma.squared.epsilon.save, type = 'l', main = paste("accept rate", out$epsilon.accept))
abline(h = s2.e)
plot(out$sigma.squared.eta.save, type = 'l', main = paste("accept rate", out$eta.accept))
abline(h = s2.s)
matplot(t(out$mu.beta.save), type = 'l')
plot(out$phi.save, type = 'l', main = paste("accept rate", out$phi.accept))
abline(h = phi)
matplot(out$fort.raster, type = 'l')
matplot(Z[, 1:100], type = 'l')
matplot((out$fort.raster - Z[, 1:100])^2, type = 'l')

