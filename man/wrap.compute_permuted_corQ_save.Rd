% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genes-meta-analyses.R
\name{wrap.compute_permuted_corQ_save}
\alias{wrap.compute_permuted_corQ_save}
\title{Wrapper to run the random shuffling of the gene correlation as several jobs on LSF}
\usage{
wrap.compute_permuted_corQ_save(DF1, DF2, method, out.df1, out.df2,
  repTimes, seed, nJobs, save = TRUE,
  importPath = "~mitra/homeshare/R321/ebits", ...)
}
\arguments{
\item{DF1}{numeric data.frame for the first study}

\item{DF2}{numeric data.frame for the second study}

\item{method}{a character string indicating which correlation coefficient is to be computed.
One of "pearson", "spearman", "kendall".}

\item{out.df1}{character string. Label for the first study. Part of the filename.}

\item{out.df2}{character string. Label for the second study. Part of the filename.}

\item{repTimes}{positive integer. Number of randomisation to create. e.g. 10,000}

\item{seed}{random seed to maintain repeatability}

\item{nJobs}{positive integer. Number of simultaneous jobs.}

\item{save}{boolean. Default: TRUE. The function tries to save in any case,
this argument defines if there is already a file with the same name
if it should be overwrite or not.}

\item{importPath}{path to ebits --- newest versions are found in \href{https://github.com/mschubert/ebits}{M. Schubert's github profile}}

\item{...}{other parameters for ebits::hpc or saveToFile}
}
\value{
Split the several randomisation of the genes correlation as different jobs on LSF
}
\description{
Wrapper to run the random shuffling of the gene correlation as several jobs on LSF
}
