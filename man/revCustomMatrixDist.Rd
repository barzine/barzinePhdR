% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-DF.R
\name{revCustomMatrixDist}
\alias{revCustomMatrixDist}
\title{Reverse a symetric data.frame (or matrix). Specifically designed for symetric distance matrix.
Diagonal is identity, i.e. the distance on the diagonal is 0
(each object is the least distant to itself)}
\usage{
revCustomMatrixDist(x, replacement = 1)
}
\arguments{
\item{x}{numeric data.frame or matrix}

\item{replacement}{numeric. Default 1 (to replace -Inf due to division by 0)}
}
\value{
numeric data.frame
}
\description{
Reverse a symetric data.frame (or matrix). Specifically designed for symetric distance matrix.
Diagonal is identity, i.e. the distance on the diagonal is 0
(each object is the least distant to itself)
}
