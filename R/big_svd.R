big_svd <- function(M,ng=1){
  p <- if(is.matrix(M)) nrow(M) else 1
  q <- if(is.matrix(M)) ncol(M) else 1
  maxnumb <- 50000
  if(p < maxnumb & q < maxnumb ){
    return(RSpectra::svds(M, nu = 1, nv = 1, k=1))
  }
  if(p > maxnumb & p/q > ng ){
    return(snm_svd(M, s = ng, h = 1))
  }
}

##-- The split and merge svd approach --##
snm_svd <- function(X, s, h){
  p <- nrow(X);  q <- ncol(X);  g <- ceiling(p/s)
  inds <- 1:g
  UMcombine<- function(UM,UMnew){
    UM$U <- c(UM$U,UMnew$U)
    UM$H <- rbind(UM$H,UMnew$H)
    return(UM)
  }
  UHmat <- foreach(i=0:(s-1), .combine = UMcombine) %dopar%{
    svdX <- svd(X[inds +g*i,], nu = min(g,q))
    D <- diag(svdX$d)
    v <- svdX$v
    list(H=D%*%t(v), U = list(svdX$u))
  }
  svdH <- RSpectra::svds(UHmat$H, h, nu =h, nv = h)
  u <- bdiag(UHmat$U)%*%svdH$u
  v <- svdH$v
  d <- svdH$d
  return(list(d=d, u=u, v=v))
}
