% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{min.lcModels}
\alias{min.lcModels}
\title{Select the lcModel with the lowest metric value}
\usage{
\method{min}{lcModels}(x, name, ...)
}
\arguments{
\item{x}{The \code{lcModels} object}

\item{name}{The name of the internal metric.}

\item{...}{Additional arguments.}
}
\value{
The lcModel with the lowest metric value
}
\description{
Select the lcModel with the lowest metric value
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
data(latrendData)
method <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time")

model1 <- latrend(method, latrendData, nClusters = 1)
model2 <- latrend(method, latrendData, nClusters = 2)
model3 <- latrend(method, latrendData, nClusters = 3)

models <- lcModels(model1, model2, model3)

min(models, "WMAE")
}
\seealso{
\link{max.lcModels} \link{externalMetric}

Other lcModels functions: 
\code{\link{as.lcModels}()},
\code{\link{lcModels}},
\code{\link{lcModels-class}},
\code{\link{max.lcModels}()},
\code{\link{plotMetric}()},
\code{\link{print.lcModels}()},
\code{\link{subset.lcModels}()}
}
\concept{lcModels functions}
