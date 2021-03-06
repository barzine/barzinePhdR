% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-prep.R
\name{maxCond.df}
\alias{maxCond.df}
\title{Give the maximum value observed for each gene
across all samples of a given tissue/condition}
\usage{
maxCond.df(DF, annot.df, cond = "all")
}
\arguments{
\item{DF}{numeric data.frame}

\item{annot.df}{annotation data.frame}

\item{cond}{vector of conditions (tissues) to be found in the returned data.frame}
}
\value{
a data.frame
}
\description{
Give the maximum value observed for each gene
across all samples of a given tissue/condition
}
