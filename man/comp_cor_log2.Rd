% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaAnalyses.R
\name{comp_cor_log2}
\alias{comp_cor_log2}
\title{Compute the correlation between two studies after their transposition
and their transformation with log2(x+1)}
\usage{
comp_cor_log2(df1, df2)
}
\arguments{
\item{df1}{numeric data.frame for the first study}

\item{df2}{numeric data.frame for the second study}
}
\value{
a named vector with the name of the first study
        and the correlation between the two datasets.
}
\description{
Compute the correlation between two studies after their transposition
and their transformation with log2(x+1)
}
