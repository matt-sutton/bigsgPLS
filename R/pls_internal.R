#' Compute Norm
#'
#' Return normed vector
#'
#' @param x matrix of numeric variables
#'
#' @return Output will be a numeric matrix or vector.
#'
my.norm <- function(x) sqrt(sum(x**2))

my.norm2 <- function(x) sum(x ^ 2)

normalize <- function(x) x / ( my.norm(x) + 1*(x==0) )

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
    y <- 1-(lambda/2)*sqrt(ji)/my.norm(vecx)
    if(y < tol) y <- 0
    res <- c(res,vecx*y)
  }
  return(res)
}

soft.thresholding.sparse.group <- function(x,ind,lambda,alpha,ind.block.zero){
  tab.ind <- c(0,ind,length(x))
  res <- NULL
  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    if(i%in%ind.block.zero) {vecx <- rep(0,ji)} else{
      temp <- soft.thresholding(vecx,lambda*alpha/2)
      vecx <- 0.5*temp*(1-lambda*(1-alpha)*sqrt(length(vecx))/sqrt(sum(temp**2)))
    }
    res <- c(res,vecx)
  }
  return(res)
}


get_lambda<- function(sparsity, ind, x, alpha){

  #-- No Sparsity --#
  if(sparsity == 0){
    return(0)
  }

  #-- sparse PLS --#
  if(is.null(ind))
    {
    lambda <- sort(abs(x))[sparsity]
    return(lambda)
  }

  #-- group PLS --#
  if( !is.null(ind) && alpha == 0 )
    {
    res <- NULL
    tab.ind <- c(0, ind, length(x))
    for (i in 1:(length(ind)+1)){
      ji <- tab.ind[i+1]-tab.ind[i]
      vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
      res <- c(res,2*my.norm(vecx)/sqrt(ji))
    }
    lambda <- sort(abs(res))[sparsity]
    return(lambda)
  }

  #-- sparse group PLS --#
  if( !is.null(ind) && alpha > 0 )
  {
    #--- ** TO DO ** ----#
  }

}
