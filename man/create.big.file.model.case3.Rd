% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_data.R
\name{create.big.file.model.case3}
\alias{create.big.file.model.case3}
\title{Simulate data matching the article}
\usage{
create.big.file.model.case3(fileX = "Xda.csv", fileY = "Yda.csv",
  size.min, p = 300, chunk.size = 9000)
}
\arguments{
\item{fileX}{filename of output X.}

\item{fileY}{filename of output Y.}

\item{size.min}{minimum size for output X.}

\item{p}{number of parameters}

\item{chunk.size}{size for chunks}
}
\value{
write an X and Y file of required size.
}
\description{
Produce a Big data file for replicaiton of the paper results.
}
