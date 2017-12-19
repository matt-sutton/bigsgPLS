Big sgPLS
=========================

`bigsgPLS` is an R package that provides an implementation of the two block PLS methods. A preliminary paper has been submitted describing the algorithm and will be made avaliable when the manuscript is accepted.

Installation
------------

Use [devtools](https://github.com/hadley/devtools) to install:

```R
library(devtools)
install_github("bigsgPLS", "matt-sutton")
```

Example Usage
-------------

```R
library(bigsgPLS)
library(bigmemory)

# Crreate some example data
# create.big.file.model.case3(size.max = 5000000000)

dataX <- read.big.matrix("Xda.csv", header = FALSE, backingfile = "Xda.bin", descriptorfile = "Xda.desc", type = "double")
Xdes <- describe(dataX)

library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

bigscale(Xdes, ng = 100)

xx <- attach.big.matrix(Xdes)
n <- nrow(xx)

dataY <- read.big.matrix("Yda.csv", header = FALSE, backingfile = "Yda.bin", descriptorfile = "Yda.desc", type = "double")

bigscale(Ydes, ng = 100)
yy <- attach.big.matrix(Ydes)

ind.block.x <- seq(100, 500, 100)
model.group.sparse.da <- algo1(Xdes, Ydes, lambda=1,regularised = "group",keepX = c(3,3),keepY = NULL,ind.block.x = ind.block.x, ind.block.y = NULL, H = 2, case = 4, epsilon = 10 ^ -6, ng = 100)


xi <- attach.big.matrix(model.group.sparse.da$xides)
omega <- attach.big.matrix(model.group.sparse.da$omegades)

which(model.group.sparse.da$loadings$X[,1]!=0)
which(model.group.sparse.da$loadings$X[,2]!=0)

xiselect <- xi[1:9000,]

par(mfrow=c(1,1))

y1 <- range(xiselect[,1])
x1 <- range(xiselect[,2])

X11(type="cairo")
par(mfrow=c(1,1),mar=c(4,4,1,1)+0.1)
plot(-4:4, -4:4, type = "n",ylim=x1,xlim=y1,xlab="Latent variable 1",ylab="Latent variable 2")
points(xiselect[1:3000,1],xiselect[1:3000,2],col="red",pch=2)
points(xiselect[3001:6000,1],xiselect[3001:6000,2],col="blue",pch=3)
points(xiselect[6001:9000,1],xiselect[6001:9000,2],col="black",pch=4)
legend("topleft",inset=0.02,c("1","2","3"),col=c("red","blue","black"),pch=c(2,3,4))

```

