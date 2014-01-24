rm(list = ls())
set.seed(1)

##
## Libraries and Subroutines
##

source('dinvgamma.R')
source('make.output.plot.R')
library(mvtnorm)
source('make.spatial.field.R')

##
## Simulate Data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 100 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.001
samp.size <- 40

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

plot.field(field$Z.list, field$H.list, locs, main = "Actual data")
plot.field(field$Y.list, field$H.list, locs)

##
## Initialize priors and tuning paramteters
##

mu.0 <- c(0, 2)#rep(0, dim(X)[2])
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
##
sigma.squared.beta.tune <- 0.025
sigma.squared.eta.tune <- 0.25
sigma.squared.epsilon.tune <- 0.075
phi.tune <- 0.45

n.mcmc <- 5000

source('mcmc.spatial.R')

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
finish <- Sys.time() - start
finish #500 iterations takes 2.23 minutes for m = 100 and reps = 100
#500 iterations takes 5.3 minutes for m = 1000 and reps = 100

##
## Plot output
##

make.output.plot(out)
## identifiability between beta_0 and sigma.squared.epsilon???
matplot(out$beta.save[1, , (n.mcmc / 10 + 1):n.mcmc], type = 'l')
matplot(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc], type = 'l')

apply(out$mu.beta.save, (n.mcmc / 10 + 1):n.mcmc], 1, mean)
dev.off()
matplot(out$fort.raster, type = 'l')


## identifiability between beta_0 and sigma.squared.epsilon???
matplot(out$beta.save[1, , (n.mcmc / 10 + 1):n.mcmc], type = 'l')
matplot(out$beta.save[2, , (n.mcmc / 10 + 1):n.mcmc], type = 'l')

coef <- matrix(nrow = reps, ncol = 2)
for(i in 1:reps){
  coef[i, ] <- lm(field$Y.list[[i]] ~ locs[field$H.list[[i]]])$coef
}
apply(coef, 2, mean)
summary(coef)
