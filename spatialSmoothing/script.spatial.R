rm(list = ls())
set.seed(1)

##
## Libraries and Subroutines
##

library(mvtnorm)
source('~/1dSpatialSim/functions/dinvgamma.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
source('~/1dSpatialSim/functions/make.spatial.field.R')
source('~/1dSpatialSim/spatialSmoothing/mcmc.spatial.R')
source('~/1dSpatialSim/functions/rMVN.R')
source('~/1dSpatialSim/plots/make.output.plot.ci.R')

# # make.output.plot <- function(out){
#   n.burn <- floor(n.mcmc / 5)
#   #x11()
#   layout(matrix(1:16, nrow = 4))
#   #
#   #   matplot(t(out$mu.beta.save)[(n.burn + 1):n.mcmc, ], type = 'l')
#   #   abline(h = beta[1], col = 'black')
#   #   abline(h = beta[2], col = 'red')
#   #
#   plot(out$sigma.squared.beta.save[(n.burn + 1):n.mcmc], type = 'l')
#   #
#   plot(out$sigma.squared.epsilon.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$epsilon.accept, 2)))
#   abline(h = s2.e, col = 'red')
#   #
#   plot(out$sigma.squared.eta.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$eta.accept, 2)))
#   abline(h = s2.s, col = 'red')
#   #
#   plot(out$phi.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$phi.accept, 2)))
#   abline(h = phi, col = 'red')
#   #
#   matplot(out$fort.raster, type = 'l')#, ylim = c(min(out$fort.raster) - 2*max(sqrt(out$var.save)), max(out$fort.raster) + 2*max(sqrt(out$var.save))))
#     matplot(out$fort.raster - 2*sqrt(out$var.save), type = 'l', add = TRUE, col = 'red', lty = 'dashed')
#     matplot(out$fort.raster + 2*sqrt(out$var.save), type = 'l', add = TRUE, col = 'red', lty = 'dashed')
#   points(X %*% beta, type = 'l', col = 'blue')  
#   #
#   plot.Z.field(field$Z.list, locs = locs, main = "True Surface")
#   #
#   #
#   plot.Y.field(field$Y.list, field$H.list, locs = locs)
#   #
#   hist(out$beta.save[1, , ][(n.burn + 1):n.mcmc])
#   abline(v = beta[1], col = 'red')
#   abline(v = quantile(out$beta.save[1, , ], probs = c(0.025, 0.975)), col = 'blue')
#   #
#   hist(out$beta.save[2, , ][(n.burn + 1):n.mcmc])
#   abline(v = beta[2], col = 'red')
#   abline(v = quantile(out$beta.save[2, , ], probs = c(0.025, 0.975)), col = 'blue')
#   #
#   MSPE <- (out$fort.raster - matrix(unlist(field$Z.list), nrow = m, byrow = FALSE))^2
#   matplot(MSPE, type = 'l', main = 'MSPE')
# }


##
## Simulate Data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 2#20 # number of spatial fields
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.01
samp.size <- 40:100

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

layout(matrix(1:2, ncol = 2))
plot.Y.field(field$Y.list[1:(reps / 2)], field$H.list[1:(reps / 2)], locs)
plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = "Full Data")

Y.list <- field$Y.list[1:(reps / 2)]
H.list <- field$H.list[1:(reps / 2)]
Z.list <- field$Z.list[(reps / 2 + 1):reps]

##
## Initialize priors and tuning paramteters
##

mu.0 <- rep(0, dim(X)[2])
sigma.squared.0 <- 25 #make sure this is large
Sigma.0 <- sigma.squared.0 * diag(dim(X)[2])
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

##
## Fit spatial MCMC kriging model
##

start <- Sys.time()
out <- mcmc.1d(field$Y.list, field$H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune)
finish <- Sys.time() - start
finish 

#5000 iterations takes 12.15 minutes for m = 1000 and reps = 10

##
## Plot output
##
#  x11()
make.output.plot(out)
MSPE.sp <- (out.sp$fort.raster - X)^2
