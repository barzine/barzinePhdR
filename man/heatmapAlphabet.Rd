% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap-wrappers.R
\name{heatmapAlphabet}
\alias{heatmapAlphabet}
\title{Draw a customised heatmap where the different variables are sorted alphabetically}
\usage{
heatmapAlphabet(data, label.A, label.B, out.A, out.B,
  col = barzinePhdR::colCbc(), digits = 1, key, margins, ...,
  report = TRUE, index = "####", annot = "diag", sorting = TRUE,
  notecol = "grey90", srtCol = 45, cexRow = 1.4, cexCol = 1.4,
  cex.axis = 1.4, mathMethod = "signif", cellAnnot)
}
\arguments{
\item{data}{correlation (numeric) data.frame to plot}

\item{label.A}{character string. Name of the first data.frame used for the correlation matrix
for the additional text.}

\item{label.B}{character string. Name of the second data.frame used for the correlation matrix
for the additional text.}

\item{out.A}{character string for the plot (x-axis).
Name of the first data.frame used for the correlation matrix}

\item{out.B}{character string for the plot (y-axis).
Name of the second data.frame used for the correlation matrix}

\item{col}{color palette to use; default barzinePhdR::colC(),
i.e. grDevices::colorRampPalette('aliceblue','darkcyan')}

\item{digits}{default: 1. How many digits should be kept for the annotation.}

\item{key}{boolean. default: TRUE. Whether the key of the heatmap should also be drawn.}

\item{margins}{numeric vector of length 2 containing the margins
(see par(mar= *)) for column and row names, respectively.}

\item{...}{other parameters that can be handled by (gplots) heatmap.2}

\item{report}{boolean. Default: TRUE. Output additional text (example for an (r)html report)}

\item{index}{character string. default:"####" Allows to customise the level of indexing through markdown tags.}

\item{annot}{character string. Chose how (if anny annotation) is added to the heatmap
'diag' (default) for the value of the diagonal;
'complete' all the values are added on the heatmap;
'custom' to use the cellAnnot given by user.
Anything else returns an heatmap that isn't annotated}

\item{sorting}{boolean; default: TRUE.
Whether the input data.frame should be rordered alphabetically before being plotted.}

\item{notecol}{character string specifying the color for cellnote text. Defaults to "grey90".}

\item{srtCol}{angle of column labels, in degrees from horizontal. Default: 45}

\item{cexRow}{positive numbers, used as cex.axis in for the row axis labeling.
Default: 1.4}

\item{cexCol}{positive numbers, used as cex.axis in for the column axis labeling.
Default: 1.4}

\item{cex.axis}{positive numbers, used as cex.axis in for the row and column axis labeling.
Default: 1.4}

\item{mathMethod}{character string that allows to pick how the correlation (in cellnote) should be rounded.
Default: "signif"}

\item{cellAnnot}{only if annot=="custom", what to be added on top of the heatmap}
}
\value{
a heatmap plot
}
\description{
Draw a customised heatmap where the different variables are sorted alphabetically
}
