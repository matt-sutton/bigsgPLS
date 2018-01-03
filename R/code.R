# create.big.file.model.case1 <- function(fileX = "X.csv", fileY = "Y.csv", size.max, p=400, q=500,chunk.size=10000) {
#   set.seed(125)
#   gam1 <- rnorm(chunk.size)
#   gam2 <- rnorm(chunk.size)
#
#   theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15),
#                 rep(0, 5), rep(-1.5, 15), rep(0, 325))
#   theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
#                 rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))
#
#   theta.y1 <- c(rep(0, 425),rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5, 15),
#                 rep(0, 5), rep(-1.5, 15))
#   theta.y2 <- c( rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
#                  rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5),rep(0, 420))
#
#   Sigmax <- matrix(0, nrow = p, ncol = p)
#   diag(Sigmax) <- sigma.e ^ 2
#   Sigmay <- matrix(0,nrow = q, ncol = q)
#   diag(Sigmay) <- sigma.e ^ 2
#
#
#   X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
#                                                                  nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, p), sigma =
#                                                                                                      Sigmax, method = "svd")
#   Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
#                                                                  nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, q), sigma =
#                                                                                                      Sigmay, method = "svd")
#   X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
#   write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
#   write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
#
#   while (file.info(fileY)$size < size.max) {
#     #print("toto")
#     gam1 <- rnorm(chunk.size)
#     gam2 <- rnorm(chunk.size)
#     X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
#                                                                    nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, p), sigma =
#                                                                                                        Sigmax, method = "svd")
#     Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
#                                                                    nrow = 2, byrow = TRUE) + rmvnorm(chunk.size, mean = rep(0, q), sigma =
#                                                                                                        Sigmay, method = "svd")
#     write.table(X, fileX, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
#     write.table(Y, fileY, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
#   }
# }
#
#
#
#
#
#
#
#
# create.big.file <- function(file = "X.csv", size.max, p, chunk.size) {
#     X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
#     write.table(X, file, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
#     while (file.info(file)$size < size.max) {
#         write.table(X, file, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
#     }
# }
#
# create.big.file2 <- function(file = "Y.csv", p, chunk.size = 1000, nb.chunk) {
#     X <- matrix(rnorm(chunk.size * p), nrow = chunk.size, ncol = p)
#     write.table(X, file, row.names = FALSE, col.names = FALSE, append = FALSE, sep = ",")
#     for (i in 1:(nb.chunk-1)) {
#         write.table(X, file, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
#     }
# }
#
# if (FALSE) {
# create.big.file(file = "X.csv", size.max=10000000,p=5,chunk.size=1000)
# dataX <- read.big.matrix("X.csv", header = FALSE, backingfile = "X.bin", descriptorfile = "X.desc", type = "double")
# Xdes <- describe(dataX)
# bigscale(Xdes, ng = 1)
# xx <- attach.big.matrix(Xdes)
# n <- nrow(xx)
#
# create.big.file2(file = "Y.csv", p=7,nb.chunk=n/1000)
# dataY <- read.big.matrix("Y.csv", header = FALSE, backingfile = "Y.bin", descriptorfile = "Y.desc", type = "double")
# Ydes <- describe(dataY)
# bigscale(Ydes, ng = 1)
# yy <- attach.big.matrix(Ydes)
# }
#
#
# if(FALSE){
# require(mixOmics)
# require(bigmemory)
# require(foreach)
# system("rm X*; rm Y*; rm xi*; rm omega*")
# data(linnerud)
# X <- linnerud$exercise
# Y <- linnerud$physiological
# linn.pls <- pls(X, Y, mode = "regression") # canonical is case2
# write.table(X,file="X.csv",col.names=F,sep=",",row.names = FALSE)
# write.table(Y,file="Y.csv",col.names=F,sep=",", row.names = FALSE)
# dataX <- read.big.matrix("X.csv", header = FALSE, backingfile = "X.bin", descriptorfile = "X.desc", type = "double")
# Xdes <- describe(dataX)
# bigscale(Xdes, ng = 1)
# xx <- attach.big.matrix(Xdes)
# n <- nrow(xx)
# dataY <- read.big.matrix("Y.csv", header = FALSE, backingfile = "Y.bin", descriptorfile = "Y.desc", type = "double")
# Ydes <- describe(dataY)
# bigscale(Ydes, ng = 1)
# yy <- attach.big.matrix(Ydes)
# res2 <- algo1(Xdes, Ydes, lambda=1, H = 3, case = 4, epsilon = 10 ^ -6, ng = 2)
# res2 <- algo1(Xdes, Ydes, sparse = TRUE,keepX = c(2,2,2),keepY = c(3,3,3),lambda=1, H = 3, case = 4, epsilon = 10 ^ -6, ng = 2)
# xi <- attach.big.matrix(res2$xides)
# # is equal to:
# linn.pls$variates$X
# omega <- attach.big.matrix(res2$omegades)
# # is equal to:
# linn.pls$variates$Y
#
# sparse <- sPLS(X,Y,ncomp = 3,keepX = c(2,2,2),mode = "regression")
# sparse2 <- spls(X,Y,ncomp = 3,keepX = c(2,2,2),mode = "regression")
# sparse$loadings$X
# sparse$loadings$Y
# sparse2$loadings$X
# sparse2$loadings$Y
#
# system("rm xi*; rm omega*")
# }
#
#
# if(FALSE){
# res4 <- algo1(Xdes, Ydes, lambda=1, H = 3, case = 2, epsilon = 10 ^ -6, ng = 2)
# xi <- attach.big.matrix(res4$xides)
# omega <- attach.big.matrix(res4$omegades)
# }
