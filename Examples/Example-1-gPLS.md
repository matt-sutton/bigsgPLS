Big sgPLS
=========================
`bigsgPLS` is an R package that provides an implementation of the two block PLS methods. A preliminary paper has been submitted describing the algorithm and will be made avaliable when the manuscript is accepted.

Big Group PLS (Test and Example)
-----------------------------------

The code below checks the method returns the same results for gPLS (using the sgPLS package).
The method is then applied to a large dataset (as used in the submitted article).

*See Also* [Example 2](Example-2-gPLS-DA.md)

```R
#--- Create some example data ---#
# First test method against the
# gPLS method from sgPLS package.
#--------------------------------#

library(mvtnorm)

set.seed(125)

n <- 100; sigma.gamma <- 1; sigma.e <- 1.5
p <- 400; q <- 500

theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15),
              rep(0, 5), rep(-1.5, 15), rep(0, 325))
theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))

theta.y1 <- c(rep(0, 425),rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5, 15),
              rep(0, 5), rep(-1.5, 15))
theta.y2 <- c( rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
               rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5),rep(0, 420))

Sigmax <- matrix(0, nrow = p, ncol = p); diag(Sigmax) <- sigma.e ^ 2
Sigmay <- matrix(0,nrow = q, ncol = q); diag(Sigmay) <- sigma.e ^ 2

set.seed(125)
gam1 <- rnorm(n); gam2 <- rnorm(n)

X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                               nrow = 2, byrow = TRUE) + rmvnorm(n, mean = rep(0, p), sigma =
                                                                                                   Sigmax, method = "svd")
Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2),
                                                               nrow = 2, byrow = TRUE) + rmvnorm(n, mean = rep(0, q), sigma =
                                                                                                   Sigmay, method = "svd")

ind.block.x <- seq(20, 380, 20)
ind.block.y <- seq(20, 480, 20)


#--------------------------------------------------------------------------#

#---- gPLS from sgPLS package ---#
library(sgPLS)

model.gPLS <- gPLS(X, Y, ncomp = 2, mode = "regression", keepX = c(4, 4),
                   keepY = c(4, 4), ind.block.x = ind.block.x , ind.block.y = ind.block.y)


result.gPLS <- select.sgpls(model.gPLS)
result.gPLS$group.size.X
result.gPLS$group.size.Y

#-- Compare to the true model
true.var.X <- c(which(theta.x1!=0),which(theta.x2!=0))
select.X <- unique(unlist(result.gPLS$select.X))
TPR <- sum(select.X%in%true.var.X)/length(true.var.X)

cat("The True Positive Rate is: ", TPR,"\n")
## The True Positive Rate is:  1

FP <- sum(select.X%in%setdiff(1:p,true.var.X))
TN <- sum(setdiff(1:p,select.X)%in%setdiff(1:p,true.var.X)) # true negatives
FPR <- FP/(FP+TN)

cat("The False Positive Rate is: ", FPR,"\n")
## The False Positive Rate is:  0.1428571

#-----------------------------------------------#

#---- gPLS from bigsgPLS package ---#

library(bigsgPLS)
library(bigmemory)

#-- Convert to bigmemory objects --#

dataX <- as.big.matrix(X)
Xdes <- describe(dataX)
bigscale(Xdes, ng = 10)

dataY <- as.big.matrix(Y)
Ydes <- describe(dataY)

bigscale(Ydes, ng = 10)

model.group.sparse <- algo1(Xdes, Ydes, regularised = "group",
                            keepX = c(4,4), keepY = c(4,4),
                            ind.block.x = ind.block.x,
                            ind.block.y = ind.block.y,
                            H = 2, case = 4, epsilon = 10 ^ -6, ng = 2)


#---- Compare the scores ---#
xi <- attach.big.matrix(model.group.sparse$xides)
omega <- attach.big.matrix(model.group.sparse$omegades)

(comparison.variatesX <- cbind(model.gPLS$variates$X,xi[1:100,]))
(comparison.variatesY <- cbind(model.gPLS$variates$Y,omega[1:100,]))

#---- Compare the weights ---#
(comparison.loadingX <- cbind(model.gPLS$loadings$X,model.group.sparse$loadings$X))
(comparison.loadingY <- cbind(model.gPLS$loadings$Y,model.group.sparse$loadings$Y))

#---- Plot to visually check against gPLS ---#
X11(type="cairo")
par(mfrow=c(2,2))
plot(model.gPLS$variates$X[,1]~xi[1:100,1],ylab="gPLS",xlab="BIG-gPLS",main="X-variates")
plot(model.gPLS$variates$Y[,1]~omega[1:100,1],ylab="sPLS",xlab="BIG-gPLS",main="Y-variates")
plot(model.gPLS$loadings$X[,1]~model.group.sparse$loadings$X[,1],ylab="gPLS",xlab="BIG-gPLS",main="X-Loading")
plot(model.gPLS$loadings$Y[,1]~model.group.sparse$loadings$Y[,1],ylab="gPLS",xlab="BIG-gPLS",main="Y-Loading")
dev.off()

#---- Plot to visually compare to truth---#
eps <- 0.2
x <- 1:p
X11(type="cairo")
par(mfrow=c(2,2),mar=c(3,3,1,1)+0.1)
norm_const <- sqrt(sum(theta.x1**2))  
plot(theta.x1~x,ylab=c(""),xlab="",col="blue",ylim=c(-1.5-eps,1.5+eps))
legend("bottomright",legend="Original value: c1",col="blue",pch=1)
plot(-norm_const*model.group.sparse$loadings$X[,1]~x,col="red",ylab=c(""),xlab="",ylim=c(-1.5-eps,1.5+eps))
legend("bottomright",legend="BIG-gPLS weight: u1",col="red",pch=1)
norm_const <- sqrt(sum(theta.y1**2))
y <- 1:q
plot(theta.y1~y,ylab=c(""),xlab="",col="blue",ylim=c(-1.5-eps,1.5+eps))
legend("bottomleft",legend="Original value: d1",col="blue",pch=1)
plot(-norm_const*model.group.sparse$loadings$Y[,1]~y,col="red",ylab=c(""),xlab="",ylim=c(-1.5-eps,1.5+eps))
legend("bottomleft",legend="BIG-gPLS weight: v1",col="red",pch=1)

#--- Method appears to return the same results ---#

#---------------------------------------------------------------#

#---- Big Data example  ---#

#--- Create some example data ---#
#
# Choose a directory for the data to be generated.
# Example from paper. Large n for X and Y matrix.
#--------------------------------#


fileX <- "../data/X.csv"
fileY <- "../data/Y.csv"

create.big.file.model.case1(fileX = fileX, fileY = fileY, size.min=5000000000, p=400, q=500,chunk.size=10000)

#-- Check the file size for the Y matrix --#
file.info(fileY)$size

library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

#-- Read the X data using bigmemory package --#
dataX <- read.big.matrix(fileX, header = FALSE, backingfile = "X.bin", descriptorfile = "X.desc", type = "double")
Xdes <- describe(dataX)

#-- Scale the X matrix --#
bigscale(Xdes, ng = 100)

#-- Read the Y data using bigmemory package --#
dataY <- read.big.matrix(fileY, header = FALSE, backingfile = "Y.bin", descriptorfile = "Y.desc", type = "double")
Ydes <- describe(dataY)

#-- Scale the Y matrix --#
bigscale(Ydes, ng = 100)

ind.block.x <- seq(20, 380, 20)
ind.block.y <- seq(20, 480, 20)

system.time(
  model.group.sparse.big <- algo1(Xdes, Ydes, regularised = "group",keepX = c(4,4),
                                  keepY = c(4,4),ind.block.x = ind.block.x ,
                                  ind.block.y = ind.block.y, H = 1, case = 4,
                                  epsilon = 10 ^ -6, ng = 100)
  )

system.time(
  model.group.sparse.big <- algo1(Xdes, Ydes, regularised = "group",keepX = c(4,4),
                                  keepY = c(4,4),ind.block.x = ind.block.x ,
                                  ind.block.y = ind.block.y, H = 2, case = 4,
                                  epsilon = 10 ^ -6, ng = 100)
)


#--- Run plots to get the weight vector estimates ---#

eps <- 0.2
x <- 1:p
X11(type="cairo")
par(mfrow=c(2,2),mar=c(3,3,1,1)+0.1)
norm_const <- sqrt(sum(theta.x1**2))  
plot(theta.x1~x,ylab=c(""),xlab="",col="blue",ylim=c(-1.5-eps,1.5+eps))
legend("bottomright",legend="Original value: c1",col="blue",pch=1)
plot(-norm_const*model.group.sparse.big$loadings$X[,1]~x,col="red",ylab=c(""),xlab="",ylim=c(-1.5-eps,1.5+eps))
legend("bottomright",legend="BIG-gPLS weight: u1",col="red",pch=1)

norm_const <- sqrt(sum(theta.y1**2))
y <- 1:q
plot(theta.y1~y,ylab=c(""),xlab="",col="blue",ylim=c(-1.5-eps,1.5+eps))
legend("bottomleft",legend="Original value: d1",col="blue",pch=1)
plot(-norm_const*model.group.sparse.big$loadings$Y[,1]~y,col="red",ylab=c(""),xlab="",ylim=c(-1.5-eps,1.5+eps))
legend("bottomleft",legend="BIG-gPLS weight: v1",col="red",pch=1)
```
