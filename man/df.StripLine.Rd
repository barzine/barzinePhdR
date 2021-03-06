% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection-filtering.R
\name{df.StripLine}
\alias{df.StripLine}
\title{Compare cell values to a threshold and fill cells
that failed the test with NA}
\usage{
df.StripLine(DF, val = 0, case = "gt")
}
\arguments{
\item{DF}{a numeric data.frame}

\item{val}{thereshold for the arithmetical comparison}

\item{case}{"ge","gt","le","lt","eq" (for df.Strip)}
}
\value{
a numeric data.frame
}
\description{
Compare cell values to a threshold and fill cells
that failed the test with NA
}
