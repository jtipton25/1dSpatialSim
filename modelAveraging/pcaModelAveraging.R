##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##

rm(list = ls())
set.seed(101)

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
X.new <- matrix(unlist(Z.list.hist), ncol = reps / 2, byrow = FALSE)
scaled <- scale.predictor(X)
# X.o <- scaled$X
X.o <- cbind(rep(1, m), scaled$X)
# X.o <- cbind(1:m, scaled$X)
# X.o <- cbind(rep(1, m), (1:m - mean(1:m)) / (sqrt(m / (m - 1)) * sd(1:m)), scaled$X)
p <- dim(X.o)[2] - 1
matplot(X, type = 'l')
matplot(X.o, type = 'l')
# X.pred <- X.new
# for(i in 1:(reps / 2)){
#   X.pred[, i] <- (X.new[, i] - scaled$scale[i, 1]) / scaled$scale[i, 2]
# }


# D <- diag(rep(max(eigen(t(X.o) %*% X.o)$values), dim(X.o)[2])) + 0.0001
# X.a <- chol(D - t(X.o) %*% X.o)

# X.c <- rbind(X.o, X.a)
# t(X.c) %*% X.c


##
## Initialize priors and tuning paramteters
##


alpha <- 2
pi.prior <- rep( 1 / 2, p)
epsilon = 0.001
n.mcmc <- 50000
lambda <- c(0, rep(1, p))
# lambda <- rep(1, p)
n.burn <- n.mcmc / 5

params <- list('vector')
params <- list(n.mcmc, alpha, pi.prior, lambda)

##
## fit mcmc using ODA model
##
# X.pca <- prcomp(X.new, center = TRUE, scale. = TRUE, retx = TRUE)$x
# 
# pca.scale <- scale.predictor(X.pca)
# X.pca.scale <- pca.scale$X
# matplot(X.pca.scale, type = 'l')

out <- mcmc.oda(Y.list = Y.list, X.o = X.o, H.list = H.list, params = params)

## Rao-blackwell estimates

beta.fit <- matrix(nrow = p + 1, ncol = reps / 2)
for(i in 1:(reps / 2)){
  for(j in 1:(p + 1)){
    if(j == 1){
#       beta.fit[j, i] <- apply(out$beta.save[j, i, ], 1, mean)
        beta.fit[j, i] <- mean(out$beta.save[j, i, ])
    } else {
#       beta.fit[j, i] <- apply(
        beta.fit[j, i] <- mean(out$rho.save[j - 1, i, ] * out$delta.save[[i]] / (out$delta.save[[i]] + lambda[2:(p + 1)][j - 1]) * out$beta.save[j - 1, i, ])
        #, 1, mean)  
    }  
  }
}

Y.pred <- matrix(nrow = m, ncol = reps / 2)
for(i in 1:(reps / 2)){
  Y.pred[, i] <- X.o %*% beta.fit[, i]  
#     Y.pred[, i] <- X.pca.scale %*% beta.fit[, i]  
}

matplot(Y.pred, type = 'l')
matplot((Y.pred - X.new)^2, type = 'l')
## mean square prediction error
MSPE <- mean((Y.pred - X.new)^2)
MSPE
# log.score <- mean(out$log.score.save[(n.burn + 1):n.mcmc])

