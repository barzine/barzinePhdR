% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-DF.R
\name{clean.infinite}
\alias{clean.infinite}
\title{Replace one value by another and remove all NA rows}
\usage{
clean.infinite(DF, motif = -Inf, token = NA)
}
\arguments{
\item{DF}{a data.frame}

\item{motif}{the value that should be changed. Default: -Inf}

\item{token}{the new value. Default NA}
}
\value{
a data.frame
}
\description{
First intention meant to change -Inf to NA
             but should fonction with other changes as well
             as long as the values to be exchanged or from equivalent classes.
}
