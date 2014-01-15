
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

plot.field <- function(Y.list, H.list, locs, main = "Observed Data", ylab = "Y", xlab = "X"){
  t <- length(Y.list)
  min.Y <- min(unlist(lapply(Y.list, min)))
  max.Y <- max(unlist(lapply(Y.list, max)))
  idx <- order(locs[H.list[[1]]])
  plot(Y.list[[1]][idx] ~ locs[H.list[[1]]][idx], type = 'l', ylim = c(min.Y, max.Y), main = main, ylab = ylab, xlab = xlab)
  if(t > 1){
    for(i in 2:t){
      idx <- order(locs[H.list[[i]]])
      lines(Y.list[[i]][idx] ~ locs[H.list[[i]]][idx], type = 'l', col = i)
    }
  }
}

##
## Initialize priors and tuning paramteters
##

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
n.mcmc <- 5000

#sigma.squared.beta.tune <- 0.25
sigma.squared.eta.tune <- 0.75
sigma.squared.epsilon.tune <- 0.25
phi.tune <- 0.75

source('mcmc.spatial.univariate.R')

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
out <- mcmc.1d(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, sigma.squared.beta, sigma.squared.eta, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta,  s.star, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
finish <- Sys.time() - start
finish #1000 iterations takes 3.75 minutes

##
## Plot output
##

#x11()
layout(matrix(1:9, nrow = 3))
plot(out$sigma.squared.beta.save, type = 'l')
plot(out$sigma.squared.epsilon.save, type = 'l', main = paste("accept rate", round(out$epsilon.accept, 2)))
abline(h = s2.e)
plot(out$sigma.squared.eta.save, type = 'l', main = paste("accept rate", round(out$eta.accept, 2)))
abline(h = s2.s)
matplot(t(out$mu.beta.save), type = 'l')
plot(out$phi.save, type = 'l', main = paste("accept rate", round(out$phi.accept, 2)))
abline(h = phi)
matplot(out$fort.raster, type = 'l')
plot.field(Z.list, H.list = list(rep(1:length(Z.list[[1]]))), locs = locs)
hist(out$mu.beta.save[1, ])
abline(v = mean(out$mu.beta.save[1, ]), col = 'red')
abline(v = quantile(out$mu.beta.save[1, ], probs = c(0.025, 0.975)), col = 'blue')
hist(out$mu.beta.save[2, ])
abline(v = mean(out$mu.beta.save[2, ]), col = 'red')
abline(v = quantile(out$mu.beta.save[2, ], probs = c(0.025, 0.975)), col = 'blue')
