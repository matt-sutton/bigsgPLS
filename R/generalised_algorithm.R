#' Unified Algorithm for sparse group PLS methods
#'
#' Function to perform big sparse group Partial Least Squares (sgPLS) in
#' the conext of datasets are divided into groups of variables.
#' The sgPLS approach enables selection at both groups and single
#' feature levels.
#'
#' @param Xdes matrix of predictors
#' @param Ydes matrix of responses
#' @param Ydes matrix of responses
#' @param lambda penalisation parameter
#' @param regularised type of regularisation
#' @param lambda penalisation parameter
#' @param keepX penalisation parameter numeric vector of length \code{H}, the number of variables
#' to keep in \eqn{X}-loadings. By default all variables are kept in the model.
#' @param keepY penalisation parameter numeric vector of length \code{H}, the number of variables
#' to keep in \eqn{Y}-loadings. By default all variables are kept in the model.
#' @param H the number of components to include in the model.
#' @param case matches the Algorithm in the paper.
#' @param ng size of chunks for fast computation
#' @param ng size of chunks for fast computation
#' @param ind.block.x a vector of integers describing the grouping of the \eqn{X}-variables.
#' @param ind.block.y a vector of integers describing the grouping of the \eqn{Y}-variables.
#' @examples
#'
#' library(bigsgPLS)
#' library(bigmemory)
#'
#' # Crreate some example data
#' # create.big.file.model.case3(size.max = 5000000000)
#'
#' dataX <- read.big.matrix("Xda.csv", header = FALSE, backingfile = "Xda.bin", descriptorfile = "Xda.desc", type = "double")
#' Xdes <- describe(dataX)
#'
#' library(doSNOW)
#' cl <- makeCluster(4)
#' registerDoSNOW(cl)
#'
#' bigscale(Xdes, ng = 100)
#'
#' xx <- attach.big.matrix(Xdes)
#' n <- nrow(xx)
#'
#' dataY <- read.big.matrix("Yda.csv", header = FALSE, backingfile = "Yda.bin", descriptorfile = "Yda.desc", type = "double")
#'
#' bigscale(Ydes, ng = 100)
#' yy <- attach.big.matrix(Ydes)
#'
#' ind.block.x <- seq(100, 500, 100)
#' model.group.sparse.da <- algo1(Xdes, Ydes, lambda=1,regularised = "group",keepX = c(3,3),keepY = NULL,ind.block.x = ind.block.x, ind.block.y = NULL, H = 2, case = 4, epsilon = 10 ^ -6, ng = 100)
#'
#' xi <- attach.big.matrix(model.group.sparse.da$xides)
#' omega <- attach.big.matrix(model.group.sparse.da$omegades)
#'
#' which(model.group.sparse.da$loadings$X[,1]!=0)
#' which(model.group.sparse.da$loadings$X[,2]!=0)
#'
#' xiselect <- xi[1:9000,]
#'
#' par(mfrow=c(1,1))
#'
#' y1 <- range(xiselect[,1])
#' x1 <- range(xiselect[,2])
#'
#' X11(type="cairo")
#' par(mfrow=c(1,1),mar=c(4,4,1,1)+0.1)
#' plot(-4:4, -4:4, type = "n",ylim=x1,xlim=y1,xlab="Latent variable 1",ylab="Latent variable 2")
#' points(xiselect[1:3000,1],xiselect[1:3000,2],col="red",pch=2)
#' points(xiselect[3001:6000,1],xiselect[3001:6000,2],col="blue",pch=3)
#' points(xiselect[6001:9000,1],xiselect[6001:9000,2],col="black",pch=4)
#' legend("topleft",inset=0.02,c("1","2","3"),col=c("red","blue","black"),pch=c(2,3,4))
#'
#'

algo1 <- function(Xdes, Ydes, regularised="none", keepX=NULL, keepY=NULL, H = 3, alpha.x = 0, alpha.y = 0,
                  case = 2, epsilon = 10 ^ -6, ng = 1, ind.block.x=NULL, ind.block.y=NULL) {

  #-- readchunk internal funciton --#
  readchunk <- function(X, g, size.chunk) {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
  }

  #-- Check data --#
  if (!(case %in% 1:4)) stop("'case' should be equal to 1, 2, 3 or 4.")

  Unew <- Vnew <- NULL
  X <- attach.big.matrix(Xdes)
  Y <- attach.big.matrix(Ydes)


  #------------------------#
  #--Set Regularisation --#

  Su <- function(v, M, lambda, ng = 1) M %*% v
  Sv <- function(u, M, lambda, ng = 1) t(M) %*% u

  if(regularised=="sparse")
    {
    #--- sparsity in terms of variables selected  ---#
    sparsity.x <- get_sparsity(keepX, ncol(X), H)
    sparsity.y <- get_sparsity(keepY, ncol(Y), H)

    #--- sPLS sparsifier ---#
    Su <- function(v, M, lambda, ng = 1) soft.thresholding(M %*% v, lambda)
    Sv <- function(u, M, lambda, ng = 1) soft.thresholding(t(M) %*% u, lambda)

  }

  if (regularised=="group")
    {

    #--- sparsity in terms of groups selected  ---#
    sparsity.x <- get_sparsity(keepX, length(ind.block.x)+1, H)
    sparsity.y <- get_sparsity(keepY, length(ind.block.y)+1, H)

    #--- gPLS sparsifier ---#
    Su <- function(v, M, lambda, ng = 1) {
      x <- M %*% v
      soft.thresholding.group(x,ind.block.x, lambda)
    }
    Sv <- function(u, M, lambda, ng = 1) {
      x <- t(M) %*% u
      soft.thresholding.group(x,ind.block.y, lambda)
    }
  }

  #--- *** sgPLS penalty under development ***  ---#

  if (regularised=="sparse group")
  {
    #--- sparsity in terms of groups selected  ---#
    sparsity.x <- get_sparsity(keepX, length(ind.block.x)+1, H)
    sparsity.y <- get_sparsity(keepY, length(ind.block.y)+1, H)

    #--- sgPLS sparsifier ---#
    Su <- function(v, M, lambda, ng = 1) {
      x <- M %*% v
      soft.thresholding.sparse.group(x, ind.block.x, lambda, alpha = alpha.x)
    }
    Sv <- function(u, M, lambda, ng = 1) {
      x <- t(M) %*% u
      soft.thresholding.sparse.group(x, ind.block.y, lambda, alpha = alpha.y)
    }
  }


  #--- set data parameters  ---#
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)

  xiH <- filebacked.big.matrix(nrow = n, ncol = H, type='double',
                               backingfile="xi.bin",
                               descriptorfile="xi.desc")

  if ((case == 2) || (case == 4)) {
    omegaH <- filebacked.big.matrix(nrow = n, ncol = H, type='double',
                                    backingfile="omega.bin",
                                    descriptorfile="omega.desc")
  }

  #-- Compute large cross product (cross product chunk)--#
  M0 <- cpc(Xdes, Ydes, ng) / (n - 1);

  if (case == 3) {
    ##Computation of A and B ## rows 2

    #-- ***** Fix this later ***** ---#

    N0y <- cpc(Ydes, Ydes, ng)
    A <- sqrtm(N0+lambda)
    B <- sqrtm(N0y+lambda)
    M0 <- A%*%M0%*%B
  }

  #-- Initialise the h - 1 vectors --#
  uhm1 <- matrix(0.0, nrow = p); vhm1 <- matrix(0.0, nrow = q)
  chm1T <- t(uhm1);  ehm1T <- t(vhm1);
  P <- Ip <- diag(1, nrow = p, ncol = p); Q <- Iq <- diag(1, nrow = q, ncol = q);

  for (h in 1:H) {
    tmp <- svd(M0, nu = 1, nv = 1)

    #-- remove sign indeterminacy --#
    i <- which.max(abs(tmp$u))
    if (tmp$u[i] <= 0) {
      tmp$u <- -tmp$u
      tmp$v <- -tmp$v
    }

    uh <- tmp$u
    vh <- tmp$v
    uprevious <- 0

    #-- Sparsifying loop for weight vectors --#
    if( (sparsity.x[h] > 0) || (sparsity.y[h] > 0) )
    {
      while ((my.norm(uh - uprevious) > epsilon)) {

        uprevious <- uh

        lambda.x <- get_lambda(sparsity.x[h], ind.block.x, M0 %*% vh, 0)
        uh <- Su(vh, M0, lambda.x, ng) ; uh <- normalize(uh)       ## row 11

        lambda.y <- get_lambda(sparsity.y[h], ind.block.y, t(M0) %*% uh, 0)
        vh <- Sv(uh, M0, lambda.y, ng) ; vh <- normalize(vh)       ## row 12
      }
    }

    #-- Compute the PLS scores --#

    xides <- describe(xiH)
    foreach(g = 1:ng) %dopar% {
      require("bigmemory")
      xiH <- attach.big.matrix(xides)
      X <- attach.big.matrix(Xdes)
      size.chunk <- nrow(xiH) / ng
      rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
      xiH[rows, h] <- X[rows,] %*% as.matrix(uh)
    }
    xih <- xiH[, h]

    omegades <- describe(omegaH)
    foreach(g = 1:ng) %dopar% {
      require("bigmemory")
      omegaH <- attach.big.matrix(omegades)
      Y <- attach.big.matrix(Ydes)
      size.chunk <- nrow(omegaH) / ng
      rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
      omegaH[rows, h] <- Y[rows,] %*% as.matrix(vh)
    }
    omegah <- omegaH[, h]


    #-- Compute the PLS adjusted weights --#

    if ( case %in% 1 ) { ## row 16
      wh <- uh ; zh <- vh ## row 17
    } ## row 18

    if ( case %in% 2 ) { ## row 19
      P <- P %*% (Ip - uhm1 %*% chm1T)  ## row 20
      Q <- Q %*% (Iq - vhm1 %*% ehm1T)  ## row 21
      wh <- P %*% uh ; zh <- Q %*% vh ## row 22
    } ## row 23

    if ( case %in% 3 ) {
      wh <- A %*% uh ; zh <- B %*% vh
      } ## row 24

    if ( case %in% 4 ) { ## row 25
      P <- P %*% (Ip - uhm1 %*% chm1T)  ## row 26
      wh <- P %*% uh  ## row 27
      zh <- vh  ## row 28
    } ## row 29


    #-- Compute the PLS loadings --#

    if ( case %in% c(1, 3) ) { ## row 30
      chm1T <- t(uh) ; ehm1T <- t(vh)  ## row 31
    } ## row 32

    if ( case %in% c(2, 4) )  chm1T <- t(xih) %*% X[,] / my.norm2(xih) ## row 33

    if ( case %in% 2 ) ehm1T <- t(omegah) %*% Y[,] / my.norm2(omegah) ## row 34

    if ( case %in% 4 ) dhm1T <- t(xih) %*% Y[,] / my.norm2(xih) ## row 35


    #-- Deflate the matrices --#

    proj <- Ip - uh %*% chm1T

    foreach(g = 1:ng, .combine = "+") %dopar% {
      X <- attach.big.matrix(Xdes)
      size.chunk <- nrow(X) / ng
      rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
      X[rows,] <- X[rows,] %*% proj
    }

    if ( case %in% 4 ) { ## row 37

      foreach(g = 1:ng) %dopar% {
        Y <- attach.big.matrix(Ydes)
        X <- attach.big.matrix(Xdes)
        size.chunk <- nrow(Y) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        Y[rows,] <- Y[rows,] - X[rows,] %*% uh %*% dhm1T
      }

      #Y <- Y[,] - xih %*% dhm1T ## row 38
      #DH[,h] <- drop(d0T) ## (row 43)

    } else { ## row 39

      proj <- Iq - vh %*% ehm1T

      foreach(g = 1:ng, .combine = "+") %dopar% {
        Y <- attach.big.matrix(Ydes)
        size.chunk <- nrow(Y) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        Y[rows,] <- Y[rows,] %*% proj
      }
      #EH[,h] <- drop(e0T) ## (row 43)

      } ## row 41

    M0 <- cpc(Xdes, Ydes, ng)

    Unew <- cbind(Unew,uh)
    Vnew <- cbind(Vnew,vh)
  }

  return(list(loadings = list(X = Unew,Y = Vnew),xides=xides,omegades=omegades,ncomp=H))
}
