% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genes-meta-analyses.R
\name{cross.sharedDistrib}
\alias{cross.sharedDistrib}
\title{Compute and compare the breath of expression of two studies and can plot them}
\usage{
cross.sharedDistrib(DF1, DF2, omit.zero = TRUE, threshold = 0,
  strict = FALSE, simplify = TRUE, Lims = c(1, 3, 10, 12),
  lapse = 3, colour, ..., annotate = TRUE, publish = TRUE,
  P1 = TRUE, P2 = FALSE, unit.DF1 = "", unit.DF2 = "",
  Prefix = "")
}
\arguments{
\item{DF1}{numeric data.frame comprising the expression data of the first study.}

\item{DF2}{numeric data.frame comprising the expression data of the second study.}

\item{omit.zero}{boolean. Default: TRUE.
 Whether the genes which do not reach the minimal expression threshold in any tissue
should be stripped from the output.}

\item{threshold}{numeric. Default: 0. Minimal level of expression to be considered as expressed.}

\item{strict}{boolean. Default:FALSE. Whether if only strict equality between the breadth of expression of each gene between the two studies should be considered
or if the "lapse" argument should be used to define 'similar' cases as well.}

\item{simplify}{boolean. Default: TRUE. Whether to focus only on the most extreme breadth of expression.
"Expressed in all" or "nearly in all" might be more descriptive than "expressed in 5" versus "expressed in 6".}

\item{Lims}{a length-4 numeric vector. Default: c(1,3,10,12).
When the analysis is in the "simplify" mode,
it allows to delimit two ranges of expression breadth to compare between the two studies.}

\item{lapse}{positive integer. Default: 3. Allow to relax constraints on perfect equality of breadth for the genes between the two studies.
A "similar" class is created. The lapse allows to define which level of similarity is acceptable.
e.g. if lapse=3, all genes that have a breadth of expression that varies at most of 3 are considered
to have similar breadth.}

\item{colour}{colour palette for the different categories.
Default: "set2" of colorbrewer2}

\item{...}{other arguments that can be used by ggplot2::theme_bw()}

\item{annotate}{boolean. Default: TRUE. Whether to annotate the figures with the counts in each category for an accurate read.}

\item{publish}{boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.}

\item{P1}{boolean. Default: TRUE. Whether to plot the figure for the first study and to return it as a port of the result.}

\item{P2}{boolean. Default: TRUE. Whether to plot the figure for the second study and to return it as a port of the result.}

\item{unit.DF1}{character string. Default: '' For the labelling, allows to input the correct unit in which the genes in the first study are expressed.}

\item{unit.DF2}{character string. Default: '' For the labelling, allows to input the correct unit in which the genes in the second study are expressed.}

\item{Prefix}{character string. Default: ''. Allows to complete the labelling on the plot.
e.g. Prefix='Expression'}
}
\value{
depending of whether P1 and P2 are true,
        the output is the data.frames with an added column for the expression breadth
        and the plot of their distribution.
}
\description{
Compute and compare the breath of expression of two studies and can plot them
}
