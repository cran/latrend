% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodLongclust.R, R/modelLongclust.R
\name{interface-longclust}
\alias{interface-longclust}
\alias{getName,lcMethodLongclust-method}
\alias{getShortName,lcMethodLongclust-method}
\alias{fit,lcMethodLongclust-method}
\alias{clusterTrajectories,lcModelLongclust-method}
\alias{postprob,lcModelLongclust-method}
\alias{converged,lcModelLongclust-method}
\alias{logLik.lcModelLongclust}
\alias{BIC.lcModelLongclust}
\title{longclust interface}
\usage{
\S4method{getName}{lcMethodLongclust}(object)

\S4method{getShortName}{lcMethodLongclust}(object)

\S4method{fit}{lcMethodLongclust}(method, data, envir, verbose, ...)

\S4method{clusterTrajectories}{lcModelLongclust}(object, at = time(object), ...)

\S4method{postprob}{lcModelLongclust}(object, ...)

\S4method{converged}{lcModelLongclust}(object, ...)

\method{logLik}{lcModelLongclust}(object, ...)

\method{BIC}{lcModelLongclust}(object, ...)
}
\arguments{
\item{object}{The object to extract the label from.}

\item{method}{The \code{lcMethod} object.}

\item{data}{The data, as a \code{data.frame}, on which the model will be trained.}

\item{envir}{The \code{environment} in which the \code{lcMethod} should be evaluated}

\item{verbose}{A \link[R.utils:Verbose]{R.utils::Verbose} object indicating the level of verbosity.}

\item{...}{Additional arguments.}

\item{at}{An optional vector, list or data frame of covariates at which to compute the cluster trajectory predictions.
If a vector is specified, this is assumed to be the time covariate. Otherwise, a named list or data frame must be provided.}
}
\description{
longclust interface
}
\seealso{
\link{lcMethodLongclust} \link[longclust]{longclust-package}
}
\keyword{internal}
