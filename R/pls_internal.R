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

soft.thresholding.sparse.group <- function(x,ind,lambda,alpha){
  tab.ind <- c(0,ind,length(x))
  res <- NULL

  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    temp <- soft.thresholding(vecx,lambda*alpha/2)

    if(my.norm(temp) <= lambda*(1-alpha)*sqrt(ji)) {
      vecx <- rep(0,ji)
      } else{

      vecx <- 0.5*temp*(1-lambda*(1-alpha)*sqrt(length(vecx))/my.norm(temp))
    }
    res <- c(res,vecx)
  }
  return(res)
}

get_sparsity<- function(keep, maxKeep, ncomp){

  #--  Compute sparsity for multiple components    -----#
  #--  Correct for mis-matched sparsity or NULL values -#

  sparsity <- maxKeep - keep

  #-- if no sparsity specified --#
  if(length(sparsity) == 0)
    {
    return(rep(0, ncomp))
  }

  #-- if sparsity mis-specified --#
  if(length(sparsity) < ncomp)
  {
    if(length(sparsity) != 1){
      print("Sparsity looks suspicious. Using a repeated sequence. \n")
    }
    return(rep(sparsity, ncomp))
  }

  return(sparsity)
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
    #-- internal function (used to set groups to 0 in sgPLS) --#
    lambda_sg_S <- function(lambda, Mv, alpha){
      return(sum(soft.thresholding(Mv,lambda*alpha/2)**2)-length(Mv)*((1-alpha)*Mv)**2)
    }

    #-- Find lambda needed for required groups --#
    res <- NULL
    tab.ind <- c(0, ind, length(x))
    for (i in 1:(length(ind)+1)){
      ji <- tab.ind[i+1]-tab.ind[i]
      vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
      upperlim <- 2*my.norm(vecx)

      res <- c(res,uniroot(lambda_sg_S, lower = 0, upper = upperlim, Mv = vecx, alpha = alpha))
    }
    lambda <- sort(res)[sparsity]
    return(lambda)
  }

}
