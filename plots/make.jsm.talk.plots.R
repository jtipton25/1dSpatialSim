##
## observed data
##

matplot(matrix(unlist(Z.list), nrow = 1000), type = 'n', main = 'Full Data', ylab = 'T', xlab = 'location')
for(i in 1:20){
  matplot(Z.list[[i]], type = 'l', col = i, add = TRUE)
}

min.Y <- min(unlist(Y.list))
max.Y <- max(unlist(Y.list))
idx <- order(locs[H.list[[1]]])
plot(Y.list[[1]][idx] ~ locs[H.list[[1]]][idx], type = 'n', ylim = c(min.Y, max.Y), main = 'Sampled data', ylab = 'T', xlab = 'location')
for(i in 1:20){
  idx <- order(locs[H.list[[i]]])
  lines(Y.list[[i]][idx] ~ locs[H.list[[i]]][idx], type = 'l', col = i)
}




## 
## pca regression
##

load('~/1dSpatialSim/data/pcaRegression.RData')

names(out)

jpeg(file = paste('~/1dSpatialSim/plots/pcaRegression', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
matplot(out$fort.raster[, 1:20], type = 'l', main = "EOF regression", ylab = "T", xlab = 'location')
dev.off()
mean((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:36])^2)

## 
## spatial regression
##

load('~/1dSpatialSim/data/spatialRegression.RData')
jpeg(file = paste('~/1dSpatialSim/plots/spatialRegression', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
matplot(out$fort.raster[, 1:20], type = 'l', main = "EOF regression with spatial random effect", ylab = "T", xlab = 'location')
dev.off()
mean((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:37])^2)
matplot((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:37])^2, type = 'l')

## 
## pp regression - 9 knots
##

load('~/1dSpatialSim/data/PPspatialRegression9knots.RData')
jpeg(file = paste('~/1dSpatialSim/plots/PPSimulation9knots', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
matplot(out$fort.raster[, 1:20], type = 'l', main = "EOF regression with predictive process - 9 knots", ylab = "T", xlab = 'location')
dev.off()
mean((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:37])^2)
  matplot((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:37])^2, type = 'l')

## 
## pp regression - 39 knots
##

load('~/1dSpatialSim/data/PPspatialRegression39knots.RData')
jpeg(file = paste('~/1dSpatialSim/plots/PPSimulation39knots', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
matplot(out$fort.raster[, 1:20], type = 'l', main = "EOF regression with predictive process - 39 knots", ylab = "T", xlab = 'location')
dev.off()
mean((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:37])^2)

## 
## pp regression - 99 knots
##

load('~/1dSpatialSim/data/PPspatialRegression39knots.RData')
jpeg(file = paste('~/1dSpatialSim/plots/PPSimulation99knots', date(), '.jpeg', sep = ''), width = 12, height = 4, quality = 100, res  = 600, units = 'in')
matplot(out$fort.raster[, 1:20], type = 'l', main = "EOF regression with predictive process - 99 knots", ylab = "T", xlab = 'location')
dev.off()
mean((out$fort.raster - matrix(unlist(field$Z.list), nrow = 1000)[, 1:37])^2)
