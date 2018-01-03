# X and Y should be standardized,  e.g.:
# Below, X should be a big.memory.matrix
readchunk <- function(X, g, size.chunk) {
  rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
  chunk <- X[rows,]
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
