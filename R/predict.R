#' predict.bigsgPLS
#'
#' Predicted values based on (sparse or sparse group) PLS models. Regression coefficient and new predictions are given using the new observations.
#'
#' @param object An object of class bigsgPLS
#' @param newX matrix or big.matrix object to make prediction on.
#' @param ng The number of chuncks used to read in the data and process using parallel computing.
#' @param comps A vector with the number of components to use in the PLS fit.
#' @param da Discriminant analysis argument to provide class estimates.
#' @param ... Further arguments passed for methods.
#'
#' @examples
#'
#' set.seed(1)
#' library(bigmemory)
#' n <- 15000
#' p <- 50
#' X = scale(matrix(rnorm(n*p), ncol = p, nrow = n))
#' y = X[,1:5] %*% 1:5 + rnorm(n)
#'
#' X.bm <- as.big.matrix(X)
#' y.bm <- as.big.matrix(y)
#'
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' getDoParWorkers()
#' fit.PLS <- bigsgpls(X.bm, y.bm, case = 4, H = 4, ng = 10, keepX = rep(5,4), regularised = "sparse")
#' pred.fit <- predict(fit.PLS, newX = X, ng = 1)
#' round(pred.fit$Beta,3)
#'
#' @export
predict.bigsgPLS <- function(object, newX, ng=1, comps = object$ncomp, da=FALSE, ...){
  x_scale <- object$scales$x_scale
  x_means <- object$scales$x_means
  y_scale <- object$scales$y_scale
  y_means <- object$scales$y_means

  comp_length <- length(comps)
  p <- ncol(newX); q <- ncol(object$Y)
  Beta <- array(0, dim = c(p,q,comp_length))
  B0 <- array(0, dim = c(p,q,comp_length))
  pred <- array(0, dim = c(nrow(newX),q,comp_length))
  for (i in 1:comp_length){
    cmps <- 1:comps[i]
    XX <- cpc(X = object$variates$X[,cmps,drop=F], Y = object$variates$X[,cmps,drop=F], ng = ng, GPU = F)
    XY <- cpc(X = object$variates$X[,cmps,drop=F], Y = object$Y, ng = ng, GPU = F)
    B <- object$adjloadings$X[,cmps,drop=F]%*%solve(XX,XY)

    B <- scale(B, center = FALSE, scale = 1 / y_scale)
    B = as.matrix(scale(t(B), center = FALSE, scale = x_scale))
    intercept = -scale(B, center = FALSE, scale = 1 / x_means)
    intercept = matrix(apply(intercept, 1, sum) + y_means, nrow = 1)
    pred[,,i] <- newX%*%t(B) + matrix(rep(1, nrow(newX)), ncol = 1) %*% intercept
    Beta[,,i] <- B
    B0[,,i] <- intercept
  }
  # prediction #
  res <- list(pred=pred, Beta = Beta, B0=intercept)
  if(da){
    res$classes = apply(pred, 1, which.max)
  }
  return(res)
}
