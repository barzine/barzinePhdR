% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genes-meta-analyses.R
\name{overlapProtTrans.Specific}
\alias{overlapProtTrans.Specific}
\title{Compute test of significance and draw venn diagram
of the observed overlap between the two data.frames.}
\usage{
overlapProtTrans.Specific(DF1, DF2, cond1, cond2, ratio = 1,
  common = TRUE, fig = TRUE, out, threshold = 1, thresholdDF1 = 0,
  categories = c("DF1", "DF2"), Col = c("coral4", "darkorchid1"),
  report = "", verbose = TRUE)
}
\arguments{
\item{DF1}{numeric data.frame of the first study to compare}

\item{DF2}{numeric data.frame of the second study to compare}

\item{cond1}{character string, column name of the first data.frame being considered for the comparison}

\item{cond2}{character string, column name of the second data.frame that being considered for the comparison}

\item{ratio}{integer. Default: 1. Possible multiplier to apply to the first data.frame genes number
to get the second data.frame genes number to consider for the comparison}

\item{common}{boolean. Default: TRUE.
Whether the two data.frames should comprise only the same genes for the comparison}

\item{fig}{boolean. Default: TRUE. Whether the figure should also be printed directly}

\item{out}{character string. "percent" or "pvalue"}

\item{threshold}{numeric. Default: 1. Minimal expression for DF2 to be considered for the comparison}

\item{thresholdDF1}{numeric. Default:0. Minimal expression for DF1 to be considered for the comparison}

\item{categories}{character string vector of two. default: c('DF1','DF2'). Allows to give the name of the studies.}

\item{Col}{numerical vectors of two. Default: "c('coral4','darkorchid1')" Colours for the venn diagram}

\item{report}{character string. Default: ''. When report='html' use cat.html instead of print}

\item{verbose}{boolean. Default: TRUE}
}
\value{
significance test and venn diagram
}
\description{
Compute test of significance and draw venn diagram
of the observed overlap between the two data.frames.
}
