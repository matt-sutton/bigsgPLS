#--     cross product chunk function    --#
#
# Return the cross product of  X : n x p   Y: n x q
# fast for large n  (p, q << n)
#
cpc <- function(X, Y, ng = 1, GPU) {

  size.chunk <- nrow(X) / ng
  if(GPU) {
    res <- foreach(g = 1:ng, .combine = "+", .packages = c("bigsgPLS","gpuR")) %dopar% {
      rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
      as.matrix(gpuR::crossprod(gpuR::gpuMatrix(X[rows,,drop=F]), gpuR::gpuMatrix(Y[rows,, drop=F]))[,])
    }
  }  else {
      res <- foreach(g = 1:ng, .combine = "+", .packages = c("bigsgPLS")) %dopar% {
        rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
        as.matrix(crossprod(X[rows,,drop=F], Y[rows,,drop=F]))
      }
  }
  return(res)
}

#--     product chunk function    --#
#
# Return the matrix product of  mat :n x p  weight: p x r
# fast for large n  (p, r << n)
#
prodchunk <- function(mat, weight, ng = 1) {
  weight <- as.matrix(weight)
  n <- nrow(mat)
  size.chunk <- n / ng

  res <- foreach(g = 1:ng, .combine = 'rbind', .packages = c("bigsgPLS")) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    mat[rows,]%*%weight
  }
  return(res)
}

#-- Internal big data function for deflating  --#
deflate <- function(X, score, loading, ng=1) {
  if(class(X) == "matrix"){
    n <- nrow(X)
    size.chunk <- n / ng
    Xnew <- foreach(g = 1:ng, .packages = c("bigsgPLS"), .combine = 'rbind') %dopar% {
      rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
      X[rows,] - score[rows]%*%loading
    }
    return(Xnew)
  } else {
    n <- nrow(X)
    size.chunk <- n / ng
    foreach(g = 1:ng, .packages = c("bigsgPLS")) %dopar% {
      rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
      X[rows,] <- X[rows,] - score[rows]%*%loading
    }
    gc()
  }
}

bigscale <- function(X, ng = 1) {

  p <- ncol(X)
  n <- nrow(X)
  size.chunk <- nrow(X) / ng
  res <- foreach(g = 1:ng, .combine = "+", .packages = c("bigsgPLS")) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,,drop=F]
    gc()
    term <- c(colSums(chunk ^ 2), colSums(chunk))
  }
  term <- res[(p + 1):(2 * p)]
  standard.deviation <- sqrt(res[1:p] / (n - 1) - term ^ 2 / (n * (n - 1)))
  average <- term / n
  remov_sd <- standard.deviation != 0

  foreach(g = 1:ng, .packages = c("bigsgPLS")) %dopar% {
    size.chunk <- nrow(X) / ng
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    X[rows,remov_sd] <- scale(X[rows,remov_sd], center = average[remov_sd], scale = standard.deviation[remov_sd])
    gc()
  }
  return(list(sd=standard.deviation, mean=average))
}


