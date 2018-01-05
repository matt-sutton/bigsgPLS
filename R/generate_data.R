#' Simulate data matching the article
#'
#' Return bigdata object for replication
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
