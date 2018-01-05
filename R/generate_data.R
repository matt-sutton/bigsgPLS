#' Simulate data matching the article
#'
#' Produce a Big data file for replicaiton of the paper results.
#'
#' @param fileX filename of output X
#' @param fileY filename of output Y
#' @param size.min minimum size for output X file
#' @param p number of parameters
#' @param chunk.size size for chunks
#' @return write an X and Y file of required size.
#'

create.big.file.model.case3 <- function(fileX = "Xda.csv", fileY = "Yda.csv", size.min, p=300,chunk.size=9000) {

  require(bigmemory)
  set.seed(123)
  mu <- c(-1,1,2)
  sigma2 <- 1
  ni <- floor(chunk.size/3)
  n <- chunk.size

  X1 <- mvtnorm::rmvnorm(ni,mean = c(rep(mu[1],p/3),rep(mu[2]+0.5,p/3),rep(0,4*p/3)),sigma = sigma2*diag(p+p))
  X2 <- mvtnorm::rmvnorm(ni,mean = c(rep(0,p/3),rep(mu[2],p/3),rep(mu[3]+0.5,p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
  X3 <- mvtnorm::rmvnorm(ni,mean = c(rep(mu[1]+0.5,p/3),rep(0,p/3),rep(mu[3],p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
  X <- rbind(X1,X2,X3)

  Y <- rep(c("1","2","3"),each=floor(n/3))
  Y = dummies::dummy(Y)

  write.table(Y,fileY,col.names=F,sep=",", append=FALSE,row.names = FALSE)
  write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")

  while (file.info(fileX)$size < size.min) {
    cat("Current file size:",file.info(fileX)$size,"\n")
    X1 <- mvtnorm::rmvnorm(ni,mean = c(rep(mu[1],p/3),rep(mu[2]+0.5,p/3),rep(0,4*p/3)),sigma = sigma2*diag(p+p))
    X2 <- mvtnorm::rmvnorm(ni,mean = c(rep(0,p/3),rep(mu[2],p/3),rep(mu[3]+0.5,p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
    X3 <- mvtnorm::rmvnorm(ni,mean = c(rep(mu[1]+0.5,p/3),rep(0,p/3),rep(mu[3],p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))

    X <- rbind(X1,X2,X3)
    Y <- rep(c("1","2","3"),each=floor(n/3))
    Y <- dummies::dummy(Y)

    write.table(Y,fileY,col.names=F,sep=",", append=TRUE,row.names = FALSE)
    write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  }
}

create.big.file.model.case1 <- function(fileX = "X.csv", fileY = "Y.csv", size.min, p=400, q=500,chunk.size=10000) {
  set.seed(125)
  sigma.gamma <- 1; sigma.e <- 1.5
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
                                                                 nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                     Sigmax, method = "svd")
  Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                 nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                     Sigmay, method = "svd")
  X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
  write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
  write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")

  while (file.info(fileY)$size < size.min) {

    cat("Current file size:",file.info(fileY)$size,"\n")
    gam1 <- rnorm(chunk.size)
    gam2 <- rnorm(chunk.size)
    X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                                   nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                       Sigmax, method = "svd")
    Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                   nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                       Sigmay, method = "svd")
    write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  }
}
