% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-prep.R
\name{extractWithinSampleAbundance}
\alias{extractWithinSampleAbundance}
\title{Extract columns that contains "WithInSampleAbundance." in their names}
\usage{
extractWithinSampleAbundance(DF, clean = TRUE)
}
\arguments{
\item{DF}{an expression data.frame}

\item{clean}{boolean; default: TRUE.
Whether rows that contains only null value should be removed}
}
\value{
a data.frame which contains only the columns that contains the tag
}
\description{
Extract columns that contains "WithInSampleAbundance." in their names
}
