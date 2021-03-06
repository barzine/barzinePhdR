% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneClassifications.R
\name{ubi.strict}
\alias{ubi.strict}
\title{Give the names of all the genes (rows) that are strictly expressed in every tissue (column) above a given threshold}
\usage{
ubi.strict(DF, threshold = 0)
}
\arguments{
\item{DF}{numeric data.frame}

\item{threshold}{numeric; default: 0. Minimal threshold above which the gene has to be observed
to be considered as expressed.}
}
\value{
List of the names of the rows complying with the threshold in all tissues (columns)
}
\description{
Give the names of all the genes (rows) that are strictly expressed in every tissue (column) above a given threshold
}
