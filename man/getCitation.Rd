% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/method.R, R/model.R
\name{getCitation}
\alias{getCitation}
\alias{getCitation,lcMethod-method}
\alias{getCitation,lcModel-method}
\title{Get citation info}
\usage{
getCitation(object, ...)

\S4method{getCitation}{lcMethod}(object, ...)

\S4method{getCitation}{lcModel}(object, ...)
}
\arguments{
\item{object}{The object}

\item{...}{Not used.}
}
\value{
A \link[utils:citation]{utils::citation} object.
}
\description{
Get a citation object indicating how to cite the underlying R packages used for estimating or representing the given method or model.
}
\seealso{
\link[utils:citation]{utils::citation}
}
