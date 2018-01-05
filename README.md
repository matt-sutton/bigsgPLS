Big sgPLS
=========================

`bigsgPLS` is an R package that provides an implementation of the two block PLS methods. A preliminary paper has been submitted describing the algorithm and will be made available when the manuscript is accepted.

Installation
------------

Use [devtools](https://github.com/hadley/devtools) to install:

```R
library(devtools)
# install_github("matt-sutton/bigsgPLS")

# -- Private Repo ---#
#
# Since the repo is private you need to set up a PAT
# in order to use the package. It's pretty easy to do.
# Just go here: https://github.com/settings/tokens
# After setting up a token, pass this into the install,
#   e.g. 

install_github("matt-sutton/bigsgPLS", auth_token = "your personal access token")

```

Example Usage
-------------

```R
library(bigsgPLS)
library(bigmemory)

#--- Create some example data ---#
#  
# Choose a directory for the data to be generated. 
# Example from paper. Large X matrix.
#--------------------------------# 

fileX <- "../data/Xda.csv"
fileY <- "../data/Yda.csv"

create.big.file.model.case3(size.min = 5000000000, chunk.size = 9000, fileX = fileX, fileY = fileY)

#-- Check the file size for the X matrix --#
file.info(fileX)$size

#-- Set up cores for parallel computation --#
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

#-- Read the X data using bigmemory package --#
dataX <- read.big.matrix(fileX, header = FALSE, backingfile = "Xda.bin", descriptorfile = "Xda.desc", type = "double")
Xdes <- describe(dataX)

#-- Scale the X matrix --#
bigscale(Xdes, ng = 100)

#-- Read the Y data using bigmemory package --#
dataY <- read.big.matrix(fileY, header = FALSE, backingfile = "Yda.bin", descriptorfile = "Yda.desc", type = "double")
Ydes <- describe(dataY)

#-- Scale the Y matrix --#
bigscale(Ydes, ng = 100)

#-- Set the block structure from the paper --#
ind.block.x <- seq(100, 500, 100)

#-- Run the Unified Algorithm with group regularisation --#
system.time(
  model.group.sparse.da <- algo1(Xdes, Ydes, regularised = "group",
                                 keepX = c(3,3),keepY = NULL,ind.block.x = ind.block.x, 
                                 ind.block.y = NULL, H = 2, case = 4, epsilon = 10 ^ -6, ng = 100)
)

# user   system  elapsed 
# 12.644    3.872   67.888 

#-- Return the PLS X and Y scores --#
xi <- attach.big.matrix(model.group.sparse.da$xides)
omega <- attach.big.matrix(model.group.sparse.da$omegades)

#-- Find the selected variables --#
which(model.group.sparse.da$loadings$X[,1]!=0)
which(model.group.sparse.da$loadings$X[,2]!=0)

#-- Get a subset of the PLS X-scores for plotting --#
xiselect <- xi[1:9000,]

par(mfrow=c(1,1))

y1 <- range(xiselect[,1])
x1 <- range(xiselect[,2])

#-- Plot the X-scores for the first two components --#

X11(type="cairo")
par(mfrow=c(1,1),mar=c(4,4,1,1)+0.1)
plot(-4:4, -4:4, type = "n",ylim=x1,xlim=y1,xlab="Latent variable 1",ylab="Latent variable 2")
points(xiselect[1:3000,1],xiselect[1:3000,2],col="red",pch=2)
points(xiselect[3001:6000,1],xiselect[3001:6000,2],col="blue",pch=3)
points(xiselect[6001:9000,1],xiselect[6001:9000,2],col="black",pch=4)
legend("topleft",inset=0.02,c("1","2","3"),col=c("red","blue","black"),pch=c(2,3,4))

```

