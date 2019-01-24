#' Unified Algorithm for sparse group PLS methods
#'
#' Function to perform big sparse group Partial Least Squares (sgPLS) in
#' the conext of datasets are divided into groups of variables.
#' The sgPLS approach enables selection at both groups and single
#' feature levels.
#'
#' @param X matrix or big.matrix object for data measured on the same samples. Corresponds to predictors in Case 2.
#' @param Y matrix or big.matrix object for data measured on the same samples. Corresponds to responses in Case 4.
#' @param regularised type of regularisation
#' @param keepX penalisation parameter numeric vector of length \code{H}, the number of variables
#' to keep in \eqn{X}-loadings. By default all variables are kept in the model.
#' @param keepY penalisation parameter numeric vector of length \code{H}, the number of variables
#' to keep in \eqn{Y}-loadings. By default all variables are kept in the model.
#' @param H the number of components to include in the model.
#' @param case matches the Algorithm in the paper.
#' @param alpha.x The mixing parameter (value between 0 and 1) related to the sparsity within group for the X dataset.
#' @param alpha.y The mixing parameter (value between 0 and 1) related to the sparsity within group for the Y dataset.
#' @param ind.block.x a vector of integers describing the grouping of the \eqn{X}-variables.
#' @param ind.block.y a vector of integers describing the grouping of the \eqn{Y}-variables.
#' @param epsilon A positive real, the tolerance used in the iterative algorithm.
#' @param ng The number of chuncks used to read in the data and process using parallel computing.
#' @param big_matrix_backing Gives the folder to use for file backed output. If NULL then output is not file backed.
#' @param GPU If TRUE then use the GPU for calculation of the chunks in the cross product. Default FALSE.
#' @param scale If TRUE then the PLS data blocks are standardized to zero means and unit variances. Default TRUE.
#' @param lambda Lambda for use in Case 3 the CCA implmenetation of PLS. Default 0.
#'
#' @examples
#'
#' set.seed(1)
#' n <- 500
#' p <- 50
#' X = scale(matrix(rnorm(n*p), ncol = p, nrow = n))
#' y = X[,1:5] %*% 1:5 + rnorm(n)
#'
#' library(bigmemory)
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
#' @import Matrix
#' @import foreach
#' @import irlba
#' @import bigmemory
#' @importFrom stats sd
#' @export

bigsgpls <- function(X,
                  Y,
                  regularised="none",
                  keepX=NULL,
                  keepY=NULL,
                  H = 3,
                  alpha.x = 0,
                  alpha.y = 0,
                  case = 2,
                  epsilon = 10 ^ -6,
                  ng = 1,
                  big_matrix_backing = NULL,
                  ind.block.x=NULL,
                  ind.block.y=NULL,
                  scale = TRUE,
                  GPU = FALSE,
                  lambda = 0) {

  if (!requireNamespace("gpuR", quietly = TRUE) && GPU) {
    stop("Package \"gpuR\" needed for the GPU computation to work. Please see https://github.com/cdeterman/gpuR/wiki for installation instructions.",
         call. = FALSE)
  }

  if(class(X) != class(Y)){
    stop("Use the same class for X and Y")
  }

  n <- nrow(X);   p <- ncol(X);   q <- ncol(Y)
  if(class(X) == "big.matrix"){
    big_matrix <- TRUE
  }
  else {
    n <- nrow(X); p <- ncol(X); q <- ncol(Y); big_matrix <- FALSE
  }

  #-- Check data --#
  if (!(case %in% 1:4)) stop("'case' should be equal to 1, 2, 3 or 4.")

  if(n > 10^5 & big_matrix == FALSE){
    cat("We recommend using bigmemory package for X and Y.")
  }

  Unew <- Vnew <- Cnew <- NULL
  Wnew <- Znew <- Enew <- NULL

  #-- Scale data if required --#
  if(!big_matrix){
    x_scale <- apply(X, 2, sd)
    x_scale[x_scale ==0] <- 1
    x_means <- colMeans(X)
    y_means <- colMeans(Y)
    y_scale <- apply(Y, 2, sd)
    y_scale[y_scale ==0] <- 1

    X <- scale(X)
    X[is.nan(X)] <- 0
    Y <- scale(Y)
  } else{
    scales <- bigscale(X, ng = ng)
    x_scale <- scales$sd
    x_scale[x_scale ==0] <- 1
    x_means <- scales$mean
    scales <- bigscale(Y, ng = ng)
    y_scale <- scales$sd
    y_scale[y_scale ==0] <- 1
    y_means <- scales$mean
  }
  Ymat <- NULL
  if(case == 4){
    Ymat <- Y[,,drop=F]
  }

  #bigscale(Y, ng = 100)

  #------------------------#
  #--Set Regularisation --#

  if(regularised=="sparse")
    {
    #--- sparsity in terms of variables selected  ---#
    sparsity.x <- get_sparsity(keepX, p, H)
    sparsity.y <- get_sparsity(keepY, q, H)

    #--- sPLS sparsifier ---#
    Su <- function(v, M, lambda) soft.thresholding(M %*% v, lambda)
    Sv <- function(u, M, lambda) soft.thresholding(t(M) %*% u, lambda)

  }

  if (regularised=="group")
    {

    #--- sparsity in terms of groups selected  ---#
    sparsity.x <- get_sparsity(keepX, length(ind.block.x)+1, H)
    sparsity.y <- get_sparsity(keepY, length(ind.block.y)+1, H)

    #--- gPLS sparsifier ---#
    Su <- function(v, M, lambda) {
      x <- M %*% v
      soft.thresholding.group(x,ind.block.x, lambda)
    }
    Sv <- function(u, M, lambda) {
      x <- t(M) %*% u
      soft.thresholding.group(x,ind.block.y, lambda)
    }
  }

  #--- sgPLS sparsifier ---#

  if (regularised=="sparse group")
  {
    #--- sparsity in terms of groups selected  ---#
    sparsity.x <- get_sparsity(keepX, length(ind.block.x)+1, H)
    sparsity.y <- get_sparsity(keepY, length(ind.block.y)+1, H)

    Su <- function(v, M, lambda) {
      x <- M %*% v
      soft.thresholding.sparse.group(x, ind.block.x, lambda, alpha = alpha.x)
    }
    Sv <- function(u, M, lambda) {
      x <- t(M) %*% u
      soft.thresholding.sparse.group(x, ind.block.y, lambda, alpha = alpha.y)
    }
  }


  #--- Initalise scores  ---#

  if(is.null(big_matrix_backing))
    {
    xiH <- matrix(nrow = n, ncol = H)
    omegaH <- matrix(nrow = n, ncol = H)

    } else {
    xiH <- bigmemory::filebacked.big.matrix(nrow = n, ncol = H, type='double',
                                 backingfile="xi.bin",
                                 descriptorfile="xi.desc", backingpath = big_matrix_backing)
    xides <- bigmemory::describe(xiH)

    omegaH <- bigmemory::filebacked.big.matrix(nrow = n, ncol = H, type='double',
                                    backingfile="omega.bin",
                                    descriptorfile="omega.desc", backingpath = big_matrix_backing)
    omegades <- bigmemory::describe(omegaH)
  }

  #-- Compute large cross product (cross product chunk)--#
  M0 <- cpc(X, Y, ng, GPU) / (n - 1);

  if (case == 3) {
    ##Computation of A and B ## rows 2

    #-- ***** Fix this later ***** ---#

    N0y <- cpc(Y, Y, ng)
    A <- expm::sqrtm(N0y+lambda)
    B <- expm::sqrtm(N0y+lambda)
    M0 <- A%*%M0%*%B
  }

  #-- Initialise the h - 1 vectors --#
  uhm1 <- matrix(0.0, nrow = p); vhm1 <- matrix(0.0, nrow = q)
  chm1T <- t(uhm1);  ehm1T <- t(vhm1);
  P <- Ip <- diag(1, nrow = p, ncol = p); Q <- Iq <- diag(1, nrow = q, ncol = q);

  for (h in 1:H) {
    tmp <- big_svd(M0, ng)

    # #-- remove sign indeterminacy --#
    i <- which.max(abs(tmp$u))
    if (tmp$u[i] <= 0) {
      tmp$u <- -tmp$u
      tmp$v <- -tmp$v
    }

    uh <- tmp$u
    vh <- tmp$v
    uprevious <- 0

    #-- Sparsifying loop for weight vectors --#
    if(regularised != "none"){
      if( (sparsity.x[h] > 0) || (sparsity.y[h] > 0) )
      {
        while ((my.norm(uh - uprevious) > epsilon)) {

          uprevious <- uh

          lambda.x <- get_lambda(sparsity.x[h], ind.block.x, M0 %*% vh, 0)
          uh <- Su(vh, M0, lambda.x) ; uh <- normalize(uh)       ## row 11

          lambda.y <- get_lambda(sparsity.y[h], ind.block.y, t(M0) %*% uh, 0)
          vh <- Sv(uh, M0, lambda.y) ; vh <- normalize(vh)       ## row 12
        }
      }
    }

    #-- Compute the PLS scores --#

    xiH[, h] <- xih <- prodchunk(X, uh, ng)
    omegaH[, h] <- omegah <- prodchunk(Y, vh, ng)

    #-- Compute the PLS adjusted weights --#

    if ( case %in% 1 ) { ## row 16
      wh <- uh ; zh <- vh ## row 17
    } ## row 18

    if ( case %in% 2 ) { ## row 19
      P <- P %*% (Ip - uh %*% chm1T)  ## row 20
      Q <- Q %*% (Iq - vh %*% ehm1T)  ## row 21
      wh <- P %*% uh ; zh <- Q %*% vh ## row 22
    } ## row 23

    if ( case %in% 3 ) {
      wh <- A %*% uh ; zh <- B %*% vh
      } ## row 24

    if ( case %in% 4 ) { ## row 25
      #P <- P %*% (Ip - uh %*% chm1T)  ## row 26
      if(h ==1) {wh=uh;zh=vh} else{
        P <- (Ip - Wnew[,1:(h-1)] %*% t(Cnew))  ## row 26
        wh <- P %*% uh  ## row 27
        zh <- vh  ## row 28
      }
    } ## row 29

    #-- Compute the PLS loadings --#

    if ( case %in% c(1, 3) ) { ## row 30
      chm1T <- t(uh) ; ehm1T <- t(vh)  ## row 31
    } ## row 32

    if ( case %in% c(2, 4) )  chm1T <- cpc( xih, X, ng, GPU) / my.norm2(xih) ## row 33

    if ( case %in% 2 ) ehm1T <- cpc(omegah, Y, ng, GPU) / my.norm2(omegah) ## row 34

    if ( case %in% 4 ) dhm1T <- cpc(xih, Y, ng, GPU) / my.norm2(xih) ## row 35

    #-- Deflate the matrices --#

    if(class(X) == "matrix") X <- deflate(X, xih, chm1T, ng = ng)  else  deflate(X, xih, chm1T, ng = ng)

    if ( case %in% 4 ) { ## row 37

      if(class(Y) == "matrix") Y <- deflate(Y, xih, dhm1T, ng = ng)  else  deflate(Y, xih, dhm1T, ng = ng)

    } else { ## row 39

      if(class(Y) == "matrix") Y <- deflate(Y, omegah, ehm1T, ng = ng)  else  deflate(Y, omegah, ehm1T, ng = ng)

      } ## row 41

    M0 <- cpc(X, Y, ng, GPU)

    Unew <- cbind(Unew,uh)
    Vnew <- cbind(Vnew,vh)
    Wnew <- cbind(Wnew,wh)
    Znew <- cbind(Znew,zh)
    Cnew <- cbind(Cnew,t(chm1T))
    Enew <- cbind(Enew,t(ehm1T))
  }

  variates <- if(!is.null(big_matrix_backing)) list(X = bigmemory::describe(xiH), Y = bigmemory::describe(omegaH)) else list(X = xiH, Y = omegaH)

  scales <- list(
    x_scale = x_scale, y_scale = y_scale, x_means = x_means, y_means = y_means
  )
  fit <- list(adjloadings = list(X = Wnew,Y = Znew),
       loadings = list(X = Unew,Y = Vnew),
       CEmat = list(Cmat = Cnew, Emat = Enew),
       variates = variates,ncomp=H,scales = scales, Y = Ymat)
  class(fit) = c(class(fit),"bigsgPLS")

  return(fit)
}
