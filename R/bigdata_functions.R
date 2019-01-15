# X and Y should be standardized,  e.g.:
# Below, X should be a big.memory.matrix
readchunk <- function(X,GPU) {
  if(GPU)
    return(gpuR::gpuMatrix(X))
  else
    return(X)
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
cpc <- function(X, Y, ng = 1, GPU) {

  X <- parse_mat(X)
  Y <- parse_mat(Y)
  if(GPU) { require(gpuR) }
  res <- foreach(g = 1:ng, .combine = "+", .packages = c("bigsgPLS")) %dopar% {
    size.chunk <- nrow(X) / ng
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk.tX <- readchunk(t(X[rows,]), GPU = GPU)
    chunk.Y <- readchunk(Y[rows,], GPU = GPU)
    as.matrix(chunk.tX[] %*% chunk.Y[])
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
  res <- foreach(g = 1:ng, .combine = 'rbind', .packages = c("bigsgPLS")) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    mat[rows,]%*%weight
  }
  return(res)
}

#-- Internal big data function for deflating  --#
deflate <- function(Xdes, score, loading, ng=1) {

  X <- parse_mat(Xdes)
  n <- nrow(X)
  size.chunk <- n / ng

  foreach(g = 1:ng, .packages = c("bigsgPLS")) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    X[rows,] <- X[rows,] - score[rows]%*%loading
    gc()
  }
  if(class(X) == "matrix")
    return(X)
}

bigscale <- function(Xdes, ng = 1) {

  X <- parse_mat(Xdes)
  p <- ncol(X)
  n <- nrow(X)

  size.chunk <- nrow(X) / ng
  res <- foreach(g = 1:ng, .combine = "+", .packages = c("bigsgPLS")) %dopar% {
    rows <- ((g - 1) * size.chunk + 1):(g * size.chunk)
    chunk <- X[rows,]
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
    X[rows,!remov_sd] <- 0
  }
  gc()
  return(list(sd=standard.deviation, mean=average))
}


