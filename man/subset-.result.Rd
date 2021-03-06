% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{[<-.result}
\alias{[<-.result}
\title{Allows to output a list directly in several objects (python like behaviour)}
\usage{
\method{[}{result}(x, ...) <- value
}
\arguments{
\item{x}{list of objects}

\item{...}{other parameters}

\item{value}{value}
}
\value{
a list of object that can be outputed directly in a list of objects
}
\description{
This function is interesting with functions
that handle several objects at the same time
The function was created by Gabor Grothendieck
who provided it graciously on the r-help list in June 2004.
https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
He originally named the function list (overriding then the built-in function).
However, I have preferred to rename it for avoiding confusion.
}
