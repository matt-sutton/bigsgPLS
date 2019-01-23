#' Simulate data matching the article
#'
#' Produce a Big data file for replicaiton of the paper results.
#'
#' @param fileX filename of output X.
#' @param fileY filename of output Y.
#' @param size.min minimum size for output X.
#' @param p number of parameters
#' @param chunk.size size for chunks
#' @return write an X and Y file of required size.
#' @export

create.big.file.model.case3 <- function(fileX = "Xda.csv", fileY = "Yda.csv", size.min, p=300,chunk.size=9000) {

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

  utils::write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
  utils::write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")

  while (file.info(fileX)$size < size.min) {
    cat("Current file size:",file.info(fileX)$size,"\n")
    X1 <- mvtnorm::rmvnorm(ni,mean = c(rep(mu[1],p/3),rep(mu[2]+0.5,p/3),rep(0,4*p/3)),sigma = sigma2*diag(p+p))
    X2 <- mvtnorm::rmvnorm(ni,mean = c(rep(0,p/3),rep(mu[2],p/3),rep(mu[3]+0.5,p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))
    X3 <- mvtnorm::rmvnorm(ni,mean = c(rep(mu[1]+0.5,p/3),rep(0,p/3),rep(mu[3],p/3),rep(0,3*p/3)),sigma = sigma2*diag(p+p))

    X <- rbind(X1,X2,X3)
    Y <- rep(c("1","2","3"),each=floor(n/3))
    Y <- dummies::dummy(Y)

    utils::write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    utils::write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  }
}

#' Simulate data matching the article
#'
#' Produce a Big data file for replicaiton of the paper results.
#'
#' @param fileX filename of output X. If NULL then no file is written and the object is returned.
#' @param fileY filename of output Y. If NULL then no file is written and the object is returned.
#' @param size.min minimum size for output X and Y given in required Gb. Default 1.
#' @param p number of parameters for X
#' @param q number of parameters for Y
#' @param chunk.size size for chunks
#' @return write an X and Y file of required size.
#' @importFrom stats rnorm
#' @export
#'
create.big.file.model.case1 <- function(fileX = NULL, fileY = NULL, size.min = 1, p=400, q=500,chunk.size=10000, ng=NULL) {
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


  memuse <- pryr::mem_change({

    X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                                 nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                     Sigmax, method = "svd");
    Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                 nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                     Sigmay, method = "svd")
  })

  numAdd <- if(is.null(ng)) ceiling(as.numeric(size.min*10^9/memuse)) - 1 else ng -1

  if(is.null(fileX) && is.null(fileY)){
    if(numAdd > 0){
      for( i in 1:numAdd ){
      X <- rbind(X,matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                                     nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                                  Sigmax, method = "svd"));
      Y <- rbind(Y,matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                     nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                                  Sigmay, method = "svd"))}
    }
    return(list(X=X, Y=Y, ind.block.x = seq(20, 380, 20), ind.block.y = seq(20, 480, 20), theta.y1=theta.y1, theta.y2=theta.y2,
                theta.x1=theta.x1, theta.x2=theta.x2))
    } else {
    utils::write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
    utils::write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
    nt <- 1
    while (file.info(fileY)$size < size.min*10^9 || numAdd < nt) {
      cat("Current file size:",file.info(fileY)$size,"\n")
      nt <- nt + 1
      if(!is.null(ng)){
        if(ng < nt) {break}
      }
      gam1 <- rnorm(chunk.size)
      gam2 <- rnorm(chunk.size)
      X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                                     nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, p), sigma =
                                                                                                                  Sigmax, method = "svd")
      Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                                     nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(chunk.size, mean = rep(0, q), sigma =
                                                                                                                  Sigmay, method = "svd")
      utils::write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
      utils::write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    }
  }
}
