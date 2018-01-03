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

algo1 <- function(Xdes, Ydes, lambda, regularised="none",keepX=NULL,keepY=NULL,H = 3,
                  case = 2, epsilon = 10 ^ -6, ng = 1,ind.block.x=NULL,ind.block.y=NULL) {

  #-- readchunk internal funciton --#
  readchunk <- function(X, g, size.chunk) {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
  }

  #-- Check data --#
  if (!(case %in% 1:4)) stop("'case' should be equal to 1, 2, 3 or 4.")
  if (!is.matrix(X)) stop("'X' should be a matrix.")
  if (!is.matrix(Y)) stop("'Y' should be a matrix.")

  Unew <- Vnew <- NULL
  X <- attach.big.matrix(Xdes)
  Y <- attach.big.matrix(Ydes)

  #------------------------#
  #--Set Regularisation --#

  Su <- function(v, M, theta.x, ng = 1) M %*% v
  Sv <- function(u, M, theta.y, ng = 1) t(M) %*% u

  if(regularised=="sparse")
    {

    sparsity.x <- ncol(X)-keepX
    sparsity.y <- ncol(Y)-keepY

    #--- sPLS  ---#
    Su <- function(v, M, theta.x, ng = 1) gsoft(M %*% v, theta.x[1])
    Sv <- function(u, M, theta.y, ng = 1) Su(u, t(M), theta.y[1], ng)

    } else if (regularised=="group")
      {
      sparsity.x <- length(ind.block.x)+1-keepX

      if(is.null(ind.block.y))
        {
        sparsity.y <- rep(0,H)
        } else {
          if (is.null(keepY)) keepY <- rep(length(ind.block.y)+1,H)
          sparsity.y <- length(ind.block.y)+1-keepY
        }
      #--- gPLS  ---#
      Su <- function(v, M, theta.x, ng = 1) {
        K <- length(theta.x) - 1
        lambda.x <- theta.x[1]
        p.vec <- theta.x[-1]
        res <- NULL
        indices <- c(0, cumsum(p.vec))
        for (k in 1:K) {
          tmp <- (crossprod.chunk(X[, (indices[k] + 1):(indices[k + 1])], Y, ng)) %*% v
          res <- cbind(res, plus.function(1.0 - 0.5 * lambda.x * sqrt(p.vec[k]) / my.norm(tmp)) * tmp)
        }
        return(res)
      }
      Sv <- function(u, M, theta.y, ng = 1) {
        return(Su(u, t(M), theta.y, ng = ng))
      }
    }

  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  xiH <- filebacked.big.matrix(nrow = n, ncol = H, type='double',
                               backingfile="xi.bin",
                               descriptorfile="xi.desc")
  if ((case == 2) || (case == 4)) {
    omegaH <- filebacked.big.matrix(nrow = n, ncol = H, type='double',
                                    backingfile="omega.bin",
                                    descriptorfile="omega.desc")
  }

  #-- compute large cross product (cross product chunk)--#
  M0 <- cpc(Xdes, Ydes, ng) / (n - 1);

  if (case == 3) {
    ##computation of A and B ## rows 2

    #-- check this with other algorithm... ---#

    N0y <- cpc(Ydes, Ydes, ng)
    A <- sqrtm(N0+lambda)
    B <- sqrtm(N0y+lambda)
    M0 <- A%*%M0%*%B
  }

  for (h in 1:H) {
    tmp <- svd(M0, nu = 1, nv = 1)

    if (remove.sign.ind) {  # so as to make the singular vectors unique (not necessary but remove sign indeterminacy)
      i <- which.max(abs(tmp$u))
      if (tmp$u[i] <= 0) {
        tmp$u <- -tmp$u
        tmp$v <- -tmp$v
      }
    }

    if (case == 3) { ### rows 14
      tmp$u <- A%*% tmp$u
      tmp$v <- B%*% tmp$uv
    }

    uold <- tmp$u
    vold <- tmp$v
    unew <- 0.0
    uprevious <- 0
    while ((my.norm(uold - uprevious) > epsilon)) {

      if(regularised=="sparse"){
        if(sparsity.x[h]==0) {lambda.x <- 0} else{
          lambda.x <- sort(abs(M0%*%matrix(vold,ncol=1)))[sparsity.x[h]]}
        unew <- soft.thresholding(M0%*%matrix(vold,ncol=1),lambda=lambda.x)
        unew <- unew/sqrt(sum(unew**2))
        ### Need to check it I put unew but  in algorithm is is uold
        if(sparsity.y[h]==0) lambda.y <- 0 else lambda.y <- sort(abs(t(M0)%*%matrix(unew,ncol=1)))[sparsity.y[h]]
        vnew <- soft.thresholding(t(M0)%*%matrix(unew,ncol=1),lambda=lambda.y)
        vnew <- vnew/sqrt(sum(vnew**2))
        uprevious <- uold
        vprevious <- vold
        uold <- unew
        vold <- vnew
      }else if (regularised=="group"){
        vecZV <- M0%*%matrix(vold,ncol=1)
        tab.ind <- c(0,ind.block.x,length(vecZV))
        res <- NULL
        for (i in 1:(length(ind.block.x)+1)){
          ji <- tab.ind[i+1]-tab.ind[i]
          vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
          res <- c(res,2*normv(vecx)/sqrt(ji))
        }
        if(sparsity.x[h]==0) lambda.x <- 0 else{
          lambda.x <- sort(res)[sparsity.x[h]]}

        unew <- soft.thresholding.group(M0%*%matrix(vold,ncol=1),ind=ind.block.x,lambda=lambda.x)
        unew <- unew/sqrt(sum(unew**2))

        if(sparsity.y[h]==0) {lambda.y <- 0} else {
          vecZV <- t(M0)%*%matrix(unew,ncol=1)
          tab.ind <- c(0,ind.block.y,length(vecZV))
          res <- NULL
          for (i in 1:(length(ind.block.y)+1)){
            ji <- tab.ind[i+1]-tab.ind[i]
            vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
            res <- c(res,2*normv(vecx)/sqrt(ji))
          }
          lambda.y <- sort(res)[sparsity.y[h]]}
        if(sparsity.y[h]==0) {vnew <- t(M0)%*%matrix(unew,ncol=1)} else {
          vnew <- soft.thresholding.group(t(M0)%*%matrix(unew,ncol=1),ind=ind.block.y,lambda=lambda.y)}
        vnew <- vnew/sqrt(sum(vnew**2))
        uprevious <- uold
        vprevious <- vold
        uold <- unew
        vold <- vnew
      }else{
        unew <- M0 %*% as.matrix(vold)
        unew <- unew / my.norm(unew)
        if (case == 3) {
          unew <-  A%*%unew
        }
        vnew <- t(M0) %*% as.matrix(uold)
        vnew <- vnew / my.norm(vnew)
        if (case == 3) {
          vnew <-  B%*%vnew
        }
        uprevious <- uold
        vprevious <- vold
        uold <- unew
        vold <- vnew
      }
    }
    if ((case ==2) || (case == 4)) {### row 26 ,
      # xih = X * unew
      xides <- describe(xiH)
      foreach(g = 1:ng) %dopar% {
        require("bigmemory")
        xiH <- attach.big.matrix(xides)
        X <- attach.big.matrix(Xdes)
        size.chunk <- nrow(xiH) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        xiH[rows, h] <- X[rows,] %*% as.matrix(unew)
      }
    } else { ### row 29 I think no need this extra else
      # xih = X * unew
      xides <- describe(xiH)
      foreach(g = 1:ng) %dopar% {
        require("bigmemory")
        xiH <- attach.big.matrix(xides)
        X <- attach.big.matrix(Xdes)
        size.chunk <- nrow(xiH) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        xiH[rows, h] <- X[rows,] %*% as.matrix(unew)
      }
    }
    if ((case == 2) || (case == 4)) { ### rows 27
      omegades <- describe(omegaH)
      foreach(g = 1:ng) %dopar% {
        require("bigmemory")
        omegaH <- attach.big.matrix(omegades)
        Y <- attach.big.matrix(Ydes)
        size.chunk <- nrow(omegaH) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        omegaH[rows, h] <- Y[rows,] %*% as.matrix(vnew)
      }
    }else{ ### rows 30
      omegades <- describe(omegaH)
      foreach(g = 1:ng) %dopar% {
        require("bigmemory")
        omegaH <- attach.big.matrix(omegades)
        Y <- attach.big.matrix(Ydes)
        size.chunk <- nrow(omegaH) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        omegaH[rows, h] <- Y[rows,] %*% as.matrix(vnew)
      }
    }

    if(h < H){

      if (case == 2) { ###row 33
        tmp <- cpc(Ydes, Ydes, ng)
        eh <- tmp %*% as.matrix(vnew) / as.vector(t(as.matrix(vnew)) %*% tmp %*% as.matrix(vnew))
      }
      if (case == 4) { ### related to rows 34
        tmpYX <- cpc(Ydes, Xdes, ng)
        tmpXX <- cpc(Xdes, Xdes, ng)
        tmpXY <- cpc(Xdes, Ydes, ng)
        uuT <- as.matrix(unew) %*% t(as.matrix(unew))
        tmp3 <- uuT %*% tmpXY
        #  dh <- tmpYX %*% as.matrix(unew) / as.vector(t(as.matrix(unew)) %*% tmpXX %*% as.matrix(unew))
        tmp4 <- as.vector(t(as.matrix(unew)) %*% tmpXX %*% as.matrix(unew))
      }

      if (case == 2) { ### rows 45
        # Yh = Yhm1 (Ip - vnew * eht)
        tmp <- diag(rep(1.0, length(vnew))) - as.matrix(vnew) %*% t(eh)
        foreach(g = 1:ng) %dopar% {
          require("bigmemory")
          Y <- attach.big.matrix(Ydes)
          size.chunk <- nrow(Y) / ng
          rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
          Y[rows,] <- Y[rows,] %*% tmp
        }
      }
      if (case == 4) { ### rows 46
        foreach(g = 1:ng) %dopar% {
          require("bigmemory")
          Y <- attach.big.matrix(Ydes)
          X <- attach.big.matrix(Xdes)
          size.chunk <- nrow(Y) / ng
          rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
          Y[rows,] <- Y[rows,] - X[rows,] %*% tmp3 / tmp4
        }

      }

      if ((case == 2) || (case == 4)) { ### rows 44
        tmp <- cpc(Xdes, Xdes, ng)
        ch <- tmp %*% as.matrix(unew) / as.vector(t(as.matrix(unew)) %*% tmp %*% as.matrix(unew))
        # Xh = Xhm1 (Ip - unew * cht)
        tmp <- diag(rep(1.0, length(unew))) - as.matrix(unew) %*% t(ch)
        foreach(g = 1:ng) %dopar% {
          require("bigmemory")
          X <- attach.big.matrix(Xdes)
          size.chunk <- nrow(X) / ng
          rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
          X[rows,] <- X[rows,] %*% tmp
        }
      }

      if ((case == 2) || (case == 4)) {
        M0 <- cpc(Xdes, Ydes, ng)
      }
    }
    Unew <- cbind(Unew,unew)
    Vnew <- cbind(Vnew,vnew)
  }

  return(list(loadings= list(X = Unew,Y = Vnew),xides=xides,omegades=omegades,ncomp=H))
}
