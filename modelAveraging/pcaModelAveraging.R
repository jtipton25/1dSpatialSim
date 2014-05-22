##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##
rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

source('~/1dSpatialSim/functions/rMVN.R')
## simulate the data
source('~/1dSpatialSim/functions/make.spatial.field.R')
## load the ODA mcmc code
source('~/1dSpatialSim/modelAveraging/mcmc.ODA.R')
## code for plotting the output
source('~/1dSpatialSim/plots/make.output.plot.ci.R')

library(statmod)
library(mvtnorm)

##
## simulate the data
##

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
X <- cbind(rep(1, m), locs)
reps <- 20 # number of spatial fields
p <- reps / 2
beta <- c(0, 2) # beta
s2.s <- 1
phi <- 0.25
s2.e <- 0.1
samp.size <- 5:40

scale.predictor <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  scale <- matrix(nrow = p, ncol = 2)
  X.tmp <- X
  for(i in 1:p){
    scale[i, ] <- c(mean(X[, i]), sqrt((n - 1) / n) * sd(X[, i]))
    X.tmp[, i] <- (X[, i] - scale[i, 1]) / scale[i, 2]
  }
  list(X = X.tmp, scale = scale)
}

field <- make.spatial.field(reps, X, beta, locs, c(s2.s, phi), method = 'exponential', s2.e, samp.size)

layout(matrix(1:2, ncol = 2))
plot.Y.field(field$Y.list[1:(reps / 2)], field$H.list[1:(reps / 2)], locs)
plot.Z.field(field$Z.list[(reps / 2 + 1):reps], locs, main = 'Full Data')

Y.list <- field$Y.list[1:(reps / 2)]
H.list <- field$H.list[1:(reps / 2)]
Z.list.hist <- field$Z.list[1:(reps / 2)]
Z.list.pca <- field$Z.list[(reps / 2 + 1):reps]
X <- matrix(unlist(Z.list.pca), ncol = reps / 2, byrow = FALSE)
scaled <- scale.predictor(X)
X.o <- scaled$X
# X.o <- cbind(rep(1, m), scaled$X)
# X.o <- cbind(1:m, scaled$X)
# X.o <- cbind(rep(1, m), 1:m, scaled$X)
matplot(X, type = 'l')
matplot(X.o, type = 'l')








D <- diag(rep(max(eigen(t(X.o) %*% X.o)$values), dim(X.o)[2])) + 0.0001
X.a <- chol(D - t(X.o) %*% X.o)
X.c <- rbind(X.o, X.a)
t(X.c) %*% X.c


##
## Initialize priors and tuning paramteters
##


alpha <- 2
pi.prior <- rep( 1 / 2, p)
epsilon = 0.001
n.mcmc <- 1000
lambda <- c(0, rep(1, p))
n.burn <- n.mcmc / 5

params <- list('vector')
params <- list(n.mcmc, alpha, pi.prior, lambda)









simdata <- make.sim.data(n.o, p, cor.vec, beta, sigma.squared, n.new)
X.o <- simdata$X.o
Y.o <- simdata$Y.o
X.new <- simdata$X.new
Y.new <- simdata$Y.new

##
## Data augmentation
##

alpha <- 2
pi.prior <- rep( 1 / 2, p)
epsilon = 0.001
n.mcmc <- 50000
lambda <- c(0, rep(1, p))
params <- list('vector')
params <- list(n.mcmc, alpha, pi.prior, lambda)
n.burn <- n.mcmc / 5
k.fold <- 8

## 
## fit mcmc using ODA model
##

out <- mcmc.oda(Y.o = Y.o, X.o = X.o, Y.new, X.new, params = params)

## Rao-blackwell estimates
beta.fit <- c(mean(out$beta.save[1, (n.burn + 1):n.mcmc]), apply(out$rho.save[, (n.burn + 1):n.mcmc] * out$delta.save / (out$delta.save + lambda[ - 1]) * out$beta.save[2:(p + 1), (n.burn + 1):n.mcmc], 1, mean))
Y.new.hat <- X.new %*% beta.fit

## mean square prediction error
MSPE <- mean((Y.new - Y.new.hat)^2)
log.score <- mean(out$log.score.save[(n.burn + 1):n.mcmc])

##
## simple linear regression with no model selection
##

mod <- lm(Y.o ~ X.o[, 2:(p + 1)])
MSPE.lm <- mean((Y.new - X.new %*% coef(mod))^2)

##
## Cross - validated selction on sigma^2_beta
## 

## make the k-fold cross-validated data
data.cv <- make.cv.data(Y.o, X.o, k.fold)

## set up parallel cluster
sfInit(parallel = TRUE, cpus = 8)
sfExportAll()
sfClusterSetupRNG()

## setup search grid for cross-validation

min.grid <- 1
max.grid <- 17
grid.size <- 16
sigma.squared.beta.cv <- seq(from = min.grid, to = max.grid, length = grid.size)

## simulate cross-validated fit

cross.validate.cv <- make.grid.search.cv(min.grid, max.grid, grid.size)
sfStop()
layout(matrix(1:2, ncol = 2))
plot(cross.validate.cv[, 2] ~ cross.validate.cv[, 1], ylab = 'MSPE', xlab = expression(sigma[beta[cv]]^2), type = 'l')
plot(cross.validate.cv[, 3] ~ cross.validate.cv[, 1], ylab = 'log score', xlab = expression(sigma[beta[cv]]^2), type = 'l')

## fit shrinkage model with cross-validated sigma.squared.beta from MSPE

sigma.squared.beta <- cross.validate.cv[, 1][which(cross.validate.cv[, 2] == min(cross.validate.cv[, 2]))]
out.cv <- mcmc.lm.cv(Y.o, X.o, Y.new, X.new, n.mcmc, sigma.squared.beta)
MSPE.cv <- mean((Y.new - apply(out.cv$y.pred.save[, (n.burn + 1):n.mcmc], 1 , mean))^2)

sigma.squared.beta <- cross.validate.cv[, 1][which(cross.validate.cv[, 3] == max(cross.validate.cv[, 3]))]
out.cv <- mcmc.lm.cv(Y.o, X.o, Y.new, X.new, n.mcmc, sigma.squared.beta)
log.score.cv <- mean(out.cv$log.score.save[(n.burn + 1):n.mcmc])

##
## fit using Lasso regression
##

## priors
alpha.epsilon <- 0.01
beta.epsilon <- 0.01
alpha.lambda <- 1
beta.lambda <- 20



out.lasso <- mcmc.lm.lasso(Y.o, X.o, Y.new, X.new, n.mcmc, alpha.epsilon, beta.epsilon, alpha.lambda, beta.lambda)
MSPE.lasso <- mean((Y.new - apply(out.lasso$y.pred.save, 1, mean))^2)
log.score.lasso <- mean(out.lasso$log.score.save[(n.burn + 1):n.mcmc])

##
## Cross-validate Lasso for prediction
##



min.grid <- 1
max.grid <- 17
grid.size <- 16

sfInit(parallel = TRUE, cpus = 8)
sfExportAll()
sfClusterSetupRNG()
sfLibrary(statmod)

cross.validate.lasso <- make.grid.search.lasso(min.grid, max.grid, grid.size)

sfStop()


##
## fit simple model with cross-validated sigma.squared.beta
##

layout(matrix(1:2, ncol = 2))
plot(cross.validate.lasso[, 2] ~ cross.validate.lasso[, 1], ylab = 'MSPE', xlab = expression(lambda[cv]^2), type = 'l')
plot(cross.validate.lasso[, 3] ~ cross.validate.lasso[, 1], ylab = 'log score', xlab = expression(lambda[cv]^2), type = 'l')

lambda.squared <- cross.validate.lasso[, 1][which(cross.validate.lasso[, 2] == min(cross.validate.lasso[, 2]))]
out.lasso.cv <- mcmc.lm.lasso.cv(Y.o, X.o, Y.new, X.new, n.mcmc, alpha.epsilon, beta.epsilon, lambda.squared)
MSPE.lasso.cv <- mean((Y.new - apply(X.new %*% out.lasso.cv$beta.save[, (n.burn + 1):n.mcmc], 1, mean))^2)

lambda.squared <- cross.validate.lasso[, 1][which(cross.validate.lasso[, 3] == max(cross.validate.lasso[, 3]))]
out.lasso.cv <- mcmc.lm.lasso.cv(Y.o, X.o, Y.new, X.new, n.mcmc, alpha.epsilon, beta.epsilon, lambda.squared)
log.score.lasso.cv <- mean(out.lasso.cv$log.score.save[(n.burn + 1):n.mcmc])

##
## save/load mcmc runs
##

# save.image('~/modelSelection/data/ODAmcmc.RData')
# load('~/modelSelection/data//ODAmcmc_May_13_2014.RData')

##
## compare MSPE from different methods
##


MSPE
MSPE.lm
MSPE.cv
MSPE.lasso
MSPE.lasso.cv

log.score
log.score.lm #need to fit a bayesian linear model
log.score.cv
log.score.lasso
log.score.lasso.cv

layout(matrix(1:2, ncol = 2))
plot(cross.validate.cv[, 2] ~ cross.validate.cv[, 1], ylab = 'MSPE', xlab = expression(sigma[beta[cv]]^2), type = 'l')
plot(cross.validate.cv[, 3] ~ cross.validate.cv[, 1], ylab = 'log score', xlab = expression(sigma[beta[cv]]^2), type = 'l')

layout(matrix(1:2, ncol = 2))
plot(cross.validate.lasso[, 2] ~ cross.validate.lasso[, 1], ylab = 'MSPE', xlab = expression(lambda[cv]^2), type = 'l')
plot(cross.validate.lasso[, 3] ~ cross.validate.lasso[, 1], ylab = 'log score', xlab = expression(lambda[cv]^2), type = 'l')

