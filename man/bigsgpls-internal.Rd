\name{bigsgPLS-internal}
\alias{cpc}
\alias{big_svd}
\alias{prodchunk}
\alias{deflate}
\title{Internal bigsgPLS functions}
\description{Internal bigsgPLS functions}
\usage{
big_svd(M,ng=1)
cpc(X,Y,ng,GPU)
prodchunk(mat,weight,ng=1)
deflate(X, score, loading, ng=1)
}
\details{
These are not intended for use by users. \code{big_svd} does svd for large datasets. \\code{cpc} Cross product using chunks for memory efficient cross products, similarly, \code{prodchunk} for matrix multiplication using chunks.  \code{deflate} deflate the matrix using chunks.
  }
\keyword{ internal }
