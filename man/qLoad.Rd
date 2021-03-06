% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{qLoad}
\alias{qLoad}
\title{Load in the working environment the object matching the name given as input
(will look for DIR.PATH in the environment)
can only load RData files}
\usage{
qLoad(objectname, path, env = globalenv())
}
\arguments{
\item{objectname}{Name of the .RData file (without the extension)
for which an object has to be created in the environment}

\item{path}{Path of the object to be loaded}

\item{env}{environment where the object should be created
Default is global environment}
}
\value{
Output directly the object inside the specified environment
}
\description{
Load in the working environment the object matching the name given as input
(will look for DIR.PATH in the environment)
can only load RData files
}
