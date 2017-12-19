#' Simulate data matching the article
#'
#' Return bigdata object for replication
#'
#' @param fileX filename of output X
#' @param fileY filename of output Y
#' @param size.max maximum size for output file
#' @param p number of parameters
#' @param chunk.size size for chunks
#' @return write an X and Y file of required size.
#'
create.big.file.model.case3 <- function(fileX = "Xda.csv", fileY = "Yda.csv", size.max, p=300,chunk.size=9000) {
  set.seed(123)
  mu <- c(-1,1,2)
  sigma2 <- 1
  ni <- chunk.size/3
  n <- chunk.size
  X1 <- rmvnorm(ni,mean = c(rep(mu[1],p/3),rep(mu[2]+0.5,p/3),rep(0,4*p/3)),sigma = sigma2*diag(p+p))
  X2 <- rmvnorm(ni,mean = c(rep(0,p/3),rep(mu[2],p/3),rep(mu[3]+0.5,p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
  X3 <- rmvnorm(ni,mean = c(rep(mu[1]+0.5,p/3),rep(0,p/3),rep(mu[3],p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
  X <- rbind(X1,X2,X3)
  Y <- rep(c("1","2","3"),each=n/3)
  YY = as.factor(Y)
  Y = unmap(as.numeric(YY))
  write.table(Y,fileY,col.names=F,sep=",", append=FALSE,row.names = FALSE)
  write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")

  while (file.info(fileX)$size < size.max) {
    X1 <- rmvnorm(ni,mean = c(rep(mu[1],p/3),rep(mu[2]+0.5,p/3),rep(0,4*p/3)),sigma = sigma2*diag(p+p))
    X2 <- rmvnorm(ni,mean = c(rep(0,p/3),rep(mu[2],p/3),rep(mu[3]+0.5,p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
    X3 <- rmvnorm(ni,mean = c(rep(mu[1]+0.5,p/3),rep(0,p/3),rep(mu[3],p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
    X <- rbind(X1,X2,X3)
    Y <- rep(c("1","2","3"),each=n/3)
    #simuData <- list(X=X, Y=factor(Y))
    #save(simuData,file="~/Dropbox/PROJECT-IEE/CODE/simuData-gPLSDA-2.RData")
    YY = as.factor(Y)
    Y = unmap(as.numeric(YY))
    write.table(Y,fileY,col.names=F,sep=",", append=TRUE,row.names = FALSE)
    write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  }
}


#' Compute Norm
#'
#' Return normed vector
#'
#' @param x matrix of numeric variables
#'
#' @return Output will be a numeric matrix or vector.
#'
normv <- function(x) sqrt(sum(x**2))


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


# X and Y should be standardized,  e.g.:
# Below, X should be a big.memory.matrix
readchunk <- function(X, g, size.chunk) {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
}

bigscale <- function(Xdes, ng = 1) {
  readchunk <- function(X, g, size.chunk) {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
  }
    res <- foreach(g = 1:ng, .combine = "+") %dopar% {
        require("bigmemory")
        X <- attach.big.matrix(Xdes)
        size.chunk <- nrow(X) / ng
        chunk <- readchunk(X, g, size.chunk)
        term <- c(colSums(chunk ^ 2), colSums(chunk))
    }
    require("bigmemory")
    X <- attach.big.matrix(Xdes)
    p <- ncol(X)
    n <- nrow(X)
    term <- res[(p + 1):(2 * p)]
    standard.deviation <- sqrt(res[1:p] / (n - 1) - term ^ 2 / (n * (n - 1)))
    average <- term / n
#    sum(x ^ 2) / (n - 1) - sum(x) ^ 2 / (n * (n - 1))
    foreach(g = 1:ng) %dopar% {
        require("bigmemory")
        X <- attach.big.matrix(Xdes)
        size.chunk <- nrow(X) / ng
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        X[rows,] <- scale(X[rows,], center = average, scale = standard.deviation)
    }
    return()
}

create.big.file.model.case1 <- function(fileX = "X.csv", fileY = "Y.csv", size.max, p=400, q=500,chunk.size=10000) {
  set.seed(125)
  gam1 <- rnorm(chunk.size)
  gam2 <- rnorm(chunk.size)

  theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15),
                rep(0, 5), rep(-1.5, 15), rep(0, 325))
  theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
                rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))

  theta.y1 <- c(rep(0, 425),rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5, 15),
                rep(0, 5), rep(-1.5, 15))
  theta.y2 <- c( rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
                 rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5),rep(0, 420))

  Sigmax <- matrix(0, nrow = p, ncol = p)
  diag(Sigmax) <- sigma.e ^ 2
  Sigmay <- matrix(0,nrow = q, ncol = q)
  diag(Sigmay) <- sigma.e ^ 2


  X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                                 nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                     Sigmax, method = "svd")
  Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                 nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                     Sigmay, method = "svd")
  X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
  write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
  write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")

  while (file.info(fileY)$size < size.max) {
    #print("toto")
    gam1 <- rnorm(chunk.size)
    gam2 <- rnorm(chunk.size)
    X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                                   nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                       Sigmax, method = "svd")
    Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                   nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                       Sigmay, method = "svd")
    write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  }
}








create.big.file <- function(file = "X.csv", size.max, p, chunk.size) {
    X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
    write.table(X, file, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
    while (file.info(file)$size < size.max) {
        write.table(X, file, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    }
}

create.big.file2 <- function(file = "Y.csv", p, chunk.size = 1000, nb.chunk) {
    X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
    write.table(X, file, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
    for (i in 1:(nb.chunk-1)) {
        write.table(X, file, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    }
}

if (FALSE) {
create.big.file(file = "X.csv", size.max=10000000,p=5,chunk.size=1000)
dataX <- read.big.matrix("X.csv", header = FALSE, backingfile = "X.bin", descriptorfile = "X.desc", type = "double")
Xdes <- describe(dataX)
bigscale(Xdes, ng = 1)
xx <- attach.big.matrix(Xdes)
n <- nrow(xx)

create.big.file2(file = "Y.csv", p=7,nb.chunk=n/1000)
dataY <- read.big.matrix("Y.csv", header = FALSE, backingfile = "Y.bin", descriptorfile = "Y.desc", type = "double")
Ydes <- describe(dataY)
bigscale(Ydes, ng = 1)
yy <- attach.big.matrix(Ydes)
}

cpc <- function(Xdes, Ydes, ng = 1) {
  readchunk <- function(X, g, size.chunk) {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
  }
    res <- foreach(g = 1:ng, .combine = "+") %dopar% {
        require("bigmemory")
        X <- attach.big.matrix(Xdes)
        Y <- attach.big.matrix(Ydes)
        size.chunk <- nrow(X) / ng
        chunk.X <- readchunk(X, g, size.chunk)
        chunk.Y <- readchunk(Y, g, size.chunk)
        term <- t(chunk.X) %*% chunk.Y
    }
    return(res)
}

my.norm <- function(x) sqrt(sum(x ^ 2))


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

algo1 <- function(Xdes, Ydes, lambda, regularised="none",keepX=NULL,keepY=NULL,H = 3, case = 2, epsilon = 10 ^ -6, ng = 1,ind.block.x=NULL,ind.block.y=NULL) {

  readchunk <- function(X, g, size.chunk) {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
  }

  Unew <- Vnew <- NULL
  X <- attach.big.matrix(Xdes)
  Y <- attach.big.matrix(Ydes)
  if(regularised=="sparse"){
    sparsity.x <- ncol(X)-keepX
    sparsity.y <- ncol(Y)-keepY
  }else if (regularised=="group"){
    sparsity.x <- length(ind.block.x)+1-keepX
    if(is.null(ind.block.y)) {sparsity.y <- rep(0,H)} else {
    if (is.null(keepY)) keepY <- rep(length(ind.block.y)+1,H)
    sparsity.y <- length(ind.block.y)+1-keepY}
  }

  n <- nrow(X)
  xiH <- filebacked.big.matrix(nrow = n, ncol = H, type='double',
                               backingfile="xi.bin",
                               descriptorfile="xi.desc")
  if ((case == 2) || (case == 4)) {
  omegaH <- filebacked.big.matrix(nrow = n, ncol = H, type='double',
                               backingfile="omega.bin",
                               descriptorfile="omega.desc")
  }
  M0 <- cpc(Xdes, Ydes, ng)
  if ((case == 3) || (case == 4)){
#    N0 <- cpc(Xdes, Xdes, ng) ## row 6,7 check it latter
  }
    if (case == 3) {
     ##computation of A and B ## rows 2
      N0y <- cpc(Ydes, Ydes, ng)
      A <- sqrtm(N0+lambda)
      B <- sqrtm(N0y+lambda)
      M0 <- A%*%M0%*%B
    }

    for (h in 1:H) {
        tmp <- svd(M0, nu = 1, nv = 1)
        i <- which.max(abs(tmp$u))
        if (tmp$u[i] <= 0) {
            tmp$u <- -tmp$u
            tmp$v <- -tmp$v
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

if(FALSE){
require(mixOmics)
require(bigmemory)
require(foreach)
system("rm X*; rm Y*; rm xi*; rm omega*")
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "regression") # canonical is case2
write.table(X,file="X.csv",col.names=F,sep=",",row.names = FALSE)
write.table(Y,file="Y.csv",col.names=F,sep=",", row.names = FALSE)
dataX <- read.big.matrix("X.csv", header = FALSE, backingfile = "X.bin", descriptorfile = "X.desc", type = "double")
Xdes <- describe(dataX)
bigscale(Xdes, ng = 1)
xx <- attach.big.matrix(Xdes)
n <- nrow(xx)
dataY <- read.big.matrix("Y.csv", header = FALSE, backingfile = "Y.bin", descriptorfile = "Y.desc", type = "double")
Ydes <- describe(dataY)
bigscale(Ydes, ng = 1)
yy <- attach.big.matrix(Ydes)
res2 <- algo1(Xdes, Ydes, lambda=1, H = 3, case = 4, epsilon = 10 ^ -6, ng = 2)
res2 <- algo1(Xdes, Ydes, sparse = TRUE,keepX = c(2,2,2),keepY = c(3,3,3),lambda=1, H = 3, case = 4, epsilon = 10 ^ -6, ng = 2)
xi <- attach.big.matrix(res2$xides)
# is equal to:
linn.pls$variates$X
omega <- attach.big.matrix(res2$omegades)
# is equal to:
linn.pls$variates$Y

sparse <- sPLS(X,Y,ncomp = 3,keepX = c(2,2,2),mode = "regression")
sparse2 <- spls(X,Y,ncomp = 3,keepX = c(2,2,2),mode = "regression")
sparse$loadings$X
sparse$loadings$Y
sparse2$loadings$X
sparse2$loadings$Y

system("rm xi*; rm omega*")
}


if(FALSE){
res4 <- algo1(Xdes, Ydes, lambda=1, H = 3, case = 2, epsilon = 10 ^ -6, ng = 2)
xi <- attach.big.matrix(res4$xides)
omega <- attach.big.matrix(res4$omegades)
}
