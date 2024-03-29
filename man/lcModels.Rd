% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{lcModels}
\alias{lcModels}
\title{Construct a list of \code{lcModel} objects}
\usage{
lcModels(...)
}
\arguments{
\item{...}{\code{lcModel}, \code{lcModels}, or a recursive \code{list} of \code{lcModel} objects. Arguments may be named.}
}
\value{
A \code{lcModels} object containing all specified \code{lcModel} objects.
}
\description{
\link[=lcModels-class]{A general overview of the lcModels class can be found here}.

The \code{lcModels()} function creates a flat (named) list of \code{lcModel} objects. Duplicates are preserved.
}
\section{Functionality}{

\itemize{
\item \link[=print.lcModels]{Print} an argument summary for each of the models.
\item \link[=as.data.frame.lcModels]{Convert} to a \code{data.frame} of method arguments.
\item \link[=subset.lcModels]{Subset} the list.
\item Compute an \link[=metric]{internal metric} or \link[=externalMetric]{external metric}.
\item Obtain the best model according to \link[=min.lcModels]{minimizing} or \link[=max.lcModels]{maximizing} a \link[=latrend-metrics]{metric}.
\item Obtain the summed \link[=estimationTime]{estimation time}.
\item \link[=plotMetric]{Plot a metric} across a variable.
\item \link[=plotClusterTrajectories]{Plot the cluster trajectories}.
}
}

\examples{

lmkmMethod <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time")
lmkmModel <- latrend(lmkmMethod, latrendData)
rngMethod <- lcMethodRandom("Y", id = "Id", time = "Time")
rngModel <- latrend(rngMethod, latrendData)

lcModels(lmkmModel, rngModel)

lcModels(defaults = c(lmkmModel, rngModel))
}
\seealso{
Other lcModels functions: 
\code{\link{as.lcModels}()},
\code{\link{lcModels-class}},
\code{\link{max.lcModels}()},
\code{\link{min.lcModels}()},
\code{\link{plotMetric}()},
\code{\link{print.lcModels}()},
\code{\link{subset.lcModels}()}
}
\concept{lcModels functions}
