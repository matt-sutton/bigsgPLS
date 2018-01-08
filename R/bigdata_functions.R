# X and Y should be standardized,  e.g.:
# Below, X should be a big.memory.matrix
readchunk <- function(X, g, size.chunk) {
  rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
  chunk <- X[rows,]
}

cpc <- function(X, Y, ng = 1) {

  res <- foreach(g = 1:ng, .combine = "+") %dopar% {
    if(class(X) == "big.matrix.descriptor") X <- bigmemory::attach.big.matrix(X)
    if(class(Y) == "big.matrix.descriptor") Y <- bigmemory::attach.big.matrix(Y)
    size.chunk <- nrow(X) / ng
    chunk.X <- readchunk(X, g, size.chunk)
    chunk.Y <- readchunk(Y, g, size.chunk)
    term <- t(chunk.X) %*% chunk.Y
  }
  return(res)
}


#-- big product --#
prodchunk <- function(des_mat, weight, ng) {

  mat <- attach.big.matrix(des_mat)
  n <- nrow(mat)
  size.chunk <- n / ng

  res <- foreach(g = 1:ng, .combine = 'rbind', .inorder = TRUE) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    mat[rows,]%*%weight
  }
  return(res)
}


bigscale <- function(Xdes, ng = 1) {

  X <- bigmemory::attach.big.matrix(Xdes)
  p <- ncol(X)
  n <- nrow(X)

  require(foreach)
  res <- foreach(g = 1:ng, .combine = "+") %dopar% {
    size.chunk <- nrow(X) / ng
    chunk <- readchunk(X, g, size.chunk)
    term <- c(colSums(chunk ^ 2), colSums(chunk))
  }
  term <- res[(p + 1):(2 * p)]
  standard.deviation <- sqrt(res[1:p] / (n - 1) - term ^ 2 / (n * (n - 1)))
  average <- term / n

  foreach(g = 1:ng) %dopar% {
    size.chunk <- nrow(X) / ng
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    X[rows,] <- scale(X[rows,], center = average, scale = standard.deviation)
  }
  return()
}
