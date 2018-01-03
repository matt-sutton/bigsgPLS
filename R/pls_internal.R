#' Compute Norm
#'
#' Return normed vector
#'
#' @param x matrix of numeric variables
#'
#' @return Output will be a numeric matrix or vector.
#'
normv <- function(x) sqrt(sum(x**2))

my.norm <- function(x) sqrt(sum(x ^ 2))

soft.thresholding <- function(x,lambda){
  tol <- .Machine$double.eps ^ 0.5
  y <- abs(x)-lambda
  test  <- y < tol
  return(sign(x)*y*(1-test))
}

soft.thresholding.group <- function(x,ind,lambda){
  tab.ind <- c(0,ind,length(x))
  tol <- .Machine$double.eps ^ 0.5
  res <- NULL
  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    y <- 1-(lambda/2)*sqrt(ji)/normv(vecx)
    if(y < tol) y <- 0
    res <- c(res,vecx*y)
  }
  return(res)
}
