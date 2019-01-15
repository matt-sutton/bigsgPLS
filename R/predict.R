predict.bigsgPLS <- function(object, newX, ng=1, comps = object$ncomp, da=FALSE){
  x_scale <- object$scales$x_scale
  x_means <- object$scales$x_means
  y_scale <- object$scales$y_scale
  y_means <- object$scales$y_means

  comp_length <- length(comps)
  p <- ncol(newX); q <- ncol(object$Y)
  Beta <- array(0, dim = c(p,q,comp_length))
  for (i in 1:comp_length){
    XX <- cpc(X = object$variates$X, Y = object$variates$X, ng = ng, GPU = F)
    XY <- cpc(X = object$variates$X, Y = object$Y, ng = ng, GPU = F)
    Beta[,,i] <- object$adjloadings$X%*%solve(XX,XY)
  }
  # prediction #
  pred <- array(0, dim = c(nrow(newX),q,comp_length))
  for( h in 1:comp_length ){
    B <- scale(Beta[,,h], center = FALSE, scale = 1 / y_scale)
    B = as.matrix(scale(t(B), center = FALSE, scale = x_scale))
    intercept = -scale(B, center = FALSE, scale = 1 / x_means)
    intercept = matrix(apply(intercept, 1, sum) + y_means, nrow = 1)
    pred[,,h] <- newX%*%t(B) + matrix(rep(1, nrow(newX)), ncol = 1) %*% intercept
  }
  res <- list(pred=pred, Beta = Beta)
  if(da){
    res$classes = apply(pred, 1, which.max)
  }
  return(res)
}
