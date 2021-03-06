% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection-filtering.R
\name{df.cutoff}
\alias{df.cutoff}
\title{Create a data.frame for which the values under a given threshold are set to 0}
\usage{
df.cutoff(DF, CUT, clean = TRUE)
}
\arguments{
\item{DF}{numeric data.frame}

\item{CUT}{The values have to be higher than CUT to be kept, otherwise value is set to 0}

\item{clean}{if TRUE, the rows where no value are above cut are filtere out}
}
\value{
a data.frame which can have a smaller number of rows if clean is TRUE
}
\description{
Create a data.frame for which the values under a given threshold are set to 0
}
