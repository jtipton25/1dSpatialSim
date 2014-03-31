##
## Plot of 1d spatial mcmc output
##

make.output.plot <- function(out){
  n.burn <- floor(n.mcmc / 10)
  #x11()
  layout(matrix(1:16, nrow = 4))
  #
  matplot(t(out$mu.beta.save)[(n.burn + 1):n.mcmc, ], type = 'l')
  abline(h = beta[1], col = 'black')
  abline(h = beta[2], col = 'red')
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
  matplot(out$fort.raster, type = 'l')
  #
  matplot(matrix(unlist(Y.pred.list), ncol = reps / 2, byrow = FALSE), type = 'l', main = 'Data to be predicted')
  #
  plot.Z.field(Z.list, locs = locs, main = "True Surface")
  #
  plot.Y.field(Y.list, H.list, locs = locs)
  #
  hist(out$mu.beta.save[1, ][(n.burn + 1):n.mcmc])
  abline(v = beta[1], col = 'red')
  abline(v = quantile(out$mu.beta.save[1, ], probs = c(0.025, 0.975)), col = 'blue')
  #
  hist(out$mu.beta.save[2, ][(n.burn + 1):n.mcmc])
  abline(v = beta[2], col = 'red')
  abline(v = quantile(out$mu.beta.save[2, ], probs = c(0.025, 0.975)), col = 'blue')
  #
  MSPE <- (out.pca$fort.raster - matrix(unlist(Y.pred.list), nrow = m, byrow = FALSE))^2
  matplot(MSPE, type = 'l', main = 'MSPE')
}

