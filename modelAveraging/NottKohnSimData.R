##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##
rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

source('~/modelSelection/functions/rMVN.R')
## simulate the data
source('~/modelSelection/functions/simNottKohn.R')
## make the cross-validation data
source('~/modelSelection/functions/make.cv.data.R')
## load the OCA mcmc code
source('~/modelSelection/ghoshAndClyde/mcmc.ODA.R')
## load shrinkage mcmc code
source('~/modelSelection/bayesianLinearRegressionCrossValidate/mcmc.lm.cv.R') # shrinkage model with shrinkage parameter chosen by cross-validation
source('~/modelSelection/functions/make.shrinkage.cv.R') # wrapper for parallelization of cross-validation of shrinkage model
## load lasso mcmc code
source('~/modelSelection/bayesianLassoRegression/fixedEffectModel/mcmc.lm.lasso.R') # lasso with shrinkage parameter in model
source('~/modelSelection/bayesianLassoRegression/fixedEffectModel/mcmc.lm.lasso.fixed.lambda.R') # lasso with shrinkage parameter fixed by cross-validation
source('~/modelSelection/functions/make.lasso.cv.R') # wrapper for paralleliztion of cross-validation of lasso model

library(snowfall)
library(rlecuyer)
library(parallel)
library(statmod)
library(mvtnorm)

##
## simulate the regression data
##


## number of covariates (not including the intercept)
p <- 15

## number of observations
n.o <- 50

## correlation vector used to correlate the first five and last five covariates\
cor.vec <- c(0.3, 0.5, 0.7, 0.9, 1.1)

## regression coefficients
beta <- c(4, 2, 0, 0, 0, -1, 0, 1.5, 0, 0, 0, 1, 0, 0.5, 0, 0)
sigma.squared <- 2.5

## number of new observations for out of sample validation
n.new <- 30

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

