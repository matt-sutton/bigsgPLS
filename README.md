Big sgPLS
=========================

`bigsgPLS` is an R package that provides an implementation of the two block PLS methods. A preliminary paper has been submitted describing the algorithm and will be made avaliable when the manuscript is accepted.

Authors
--------
This is joint work with Pierre Lafeye De Micheaux and Benoit Liquet. A preliminary report describing the PLS methods and some of the statistical properties is avaliable on [ArXiv Pre-prints.](https://arxiv.org/abs/1702.07066)

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

Examples from the bigsgPLS paper can be replcated by following the process described in: [Example 1](Examples/Example-1-gPLS.md) and [Example 2](Examples/Example-2-gPLS-DA.md)
