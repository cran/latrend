% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelPartition.R, R/methodFunction.R,
%   R/methodRandom.R, R/modelStratify.R, R/modelWeightedPartition.R
\name{clusterTrajectories,lcModelPartition-method}
\alias{clusterTrajectories,lcModelPartition-method}
\alias{converged,lcModelPartition-method}
\alias{getName,lcModelPartition-method}
\alias{getShortName,lcModelPartition-method}
\alias{postprob,lcModelPartition-method}
\alias{interface-custom}
\alias{getArgumentDefaults,lcMethodFunction-method}
\alias{getName,lcMethodFunction-method}
\alias{getShortName,lcMethodFunction-method}
\alias{prepareData,lcMethodFunction-method}
\alias{fit,lcMethodFunction-method}
\alias{getArgumentDefaults,lcMethodRandom-method}
\alias{getName,lcMethodRandom-method}
\alias{getShortName,lcMethodRandom-method}
\alias{fit,lcMethodRandom-method}
\alias{converged,lcModelStratify-method}
\alias{predictPostprob,lcModelStratify-method}
\alias{clusterTrajectories,lcModelWeightedPartition-method}
\alias{converged,lcModelWeightedPartition-method}
\alias{getName,lcModelWeightedPartition-method}
\alias{getShortName,lcModelWeightedPartition-method}
\alias{postprob,lcModelWeightedPartition-method}
\title{function interface}
\usage{
\S4method{clusterTrajectories}{lcModelPartition}(
  object,
  at = time(object),
  center = object@center,
  approxFun = approx,
  ...
)

\S4method{converged}{lcModelPartition}(object, ...)

\S4method{getName}{lcModelPartition}(object, ...)

\S4method{getShortName}{lcModelPartition}(object, ...)

\S4method{postprob}{lcModelPartition}(object, ...)

\S4method{getArgumentDefaults}{lcMethodFunction}(object)

\S4method{getName}{lcMethodFunction}(object)

\S4method{getShortName}{lcMethodFunction}(object)

\S4method{prepareData}{lcMethodFunction}(method, data, verbose)

\S4method{fit}{lcMethodFunction}(method, data, envir, verbose)

\S4method{getArgumentDefaults}{lcMethodRandom}(object)

\S4method{getName}{lcMethodRandom}(object)

\S4method{getShortName}{lcMethodRandom}(object)

\S4method{fit}{lcMethodRandom}(method, data, envir, verbose, ...)

\S4method{converged}{lcModelStratify}(object, ...)

\S4method{predictPostprob}{lcModelStratify}(object, newdata = NULL, ...)

\S4method{clusterTrajectories}{lcModelWeightedPartition}(
  object,
  at = time(object),
  center = weighted.meanNA,
  approxFun = approx,
  ...
)

\S4method{converged}{lcModelWeightedPartition}(object, ...)

\S4method{getName}{lcModelWeightedPartition}(object, ...)

\S4method{getShortName}{lcModelWeightedPartition}(object, ...)

\S4method{postprob}{lcModelWeightedPartition}(object, ...)
}
\arguments{
\item{object}{The object.}

\item{center}{The function to use to compute the cluster trajectory center at the respective moment in time.}

\item{...}{Not used.}

\item{method}{An object inheriting from \code{lcMethod} with all its arguments having been evaluated and finalized.}

\item{data}{A \code{data.frame} representing the transformed training data.}

\item{verbose}{A \link[R.utils:Verbose]{R.utils::Verbose} object indicating the level of verbosity.}

\item{envir}{The \code{environment} containing variables generated by \code{\link[=prepareData]{prepareData()}} and \code{\link[=preFit]{preFit()}}.}

\item{newdata}{Optional \code{data.frame} for which to compute the posterior probability. If omitted, the model training data is used.}
}
\description{
function interface
}
\seealso{
\link{lcMethodRandom} \link{lcMethodStratify} \link{lcModelPartition} \link{lcModelWeightedPartition}
}
\keyword{internal}
