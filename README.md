Big sgPLS
=========================

`bigsgPLS` is an R package that provides an implementation of the two block PLS methods. The method makes use of bigmemory and matrix algebra by chunks to deal with datasets too large for R.

Authors
--------
This is joint work with Pierre Lafeye De Micheaux and Benoit Liquet. A preliminary paper describing the PLS methods and some of the statistical properties is avaliable on [ArXiv Pre-prints.](https://arxiv.org/abs/1702.07066)

Installation
------------

Use [devtools](https://github.com/hadley/devtools) to install:

```R
library(devtools)
install_github("matt-sutton/bigsgPLS", host = "https://api.github.com")
```

Example Usage
-------------

Examples from the bigsgPLS paper can be replicated by following the process described in: [Example 1](Examples/Example-1-gPLS.md) and [Example 2](Examples/Example-2-gPLS-DA.md) and [Example 3](Examples/Example-3-PLS.md)
