# X and Y should be standardized,  e.g.:
# Below, X should be a big.memory.matrix
readchunk <- function(X, g, size.chunk) {
  rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
  chunk <- X[rows,]
}

#-- Parse either matrix or descriptor --#
parse_mat <- function(X){

  X <- switch(class(X),
         "big.matrix.descriptor" = bigmemory::attach.big.matrix(X),
         "matrix" = X,
         stop("data should be a matrix or big matrix descriptor.")
  )
}


#--     cross product chunk function    --#
#
# Return the cross product of  X : n x p   Y: n x q
# fast for large n  (p, q << n)
#
cpc <- function(X, Y, ng = 1) {

  X <- parse_mat(X)
  Y <- parse_mat(Y)
  require(foreach)
  res <- foreach(g = 1:ng, .combine = "+") %dopar% {
    size.chunk <- nrow(X) / ng
    chunk.X <- readchunk(X, g, size.chunk)
    chunk.Y <- readchunk(Y, g, size.chunk)
    term <- t(chunk.X) %*% chunk.Y
  }
  return(res)
}


#--     product chunk function    --#
#
# Return the matrix product of  mat :n x p  weight: p x r
# fast for large n  (p, r << n)
#
prodchunk <- function(des_mat, weight, ng) {

  mat <- parse_mat(des_mat)
  n <- nrow(mat)
  size.chunk <- n / ng

  require(foreach)
  res <- foreach(g = 1:ng, .combine = 'rbind', .inorder = TRUE) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    mat[rows,]%*%weight
  }
  return(res)
}


bigscale <- function(Xdes, ng = 1) {

  X <- parse_mat(Xdes)
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
