% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-transform.R
\name{transformPredict}
\alias{transformPredict}
\alias{transformPredict,NULL,lcModel-method}
\alias{transformPredict,vector,lcModel-method}
\alias{transformPredict,matrix,lcModel-method}
\alias{transformPredict,data.frame,lcModel-method}
\title{Helper function for custom lcModel classes implementing predict.lcModel()}
\usage{
transformPredict(pred, model, newdata)

\S4method{transformPredict}{NULL,lcModel}(pred, model, newdata)

\S4method{transformPredict}{vector,lcModel}(pred, model, newdata)

\S4method{transformPredict}{matrix,lcModel}(pred, model, newdata)

\S4method{transformPredict}{data.frame,lcModel}(pred, model, newdata)
}
\arguments{
\item{pred}{The (per-cluster) predictions for \code{newdata}.}

\item{model}{The \code{lcModel} for which the prediction was made.}

\item{newdata}{A \code{data.frame} containing the input data to predict for.}
}
\value{
A \code{data.frame} with the predictions, or a list of cluster-specific prediction \code{data.frame}s.
}
\description{
A helper function for implementing the \link[=predict.lcModel]{predict.lcModel()} method as part of your own \code{lcModel} class, ensuring the correct output type and format (see the Value section).
Note that this function has no use outside of ensuring valid output for \code{predict.lcModel}.
For implementing \code{lcModel} predictions from scratch, it is advisable to implement \link{predictForCluster} instead of \link{predict.lcModel}.

The prediction ordering corresponds to the observation ordering of the \code{newdata} argument.

By default, \code{transformPredict()} accepts one of the following inputs:
\describe{
\item{\code{data.frame}}{A \code{data.frame} in long format providing a cluster-specific prediction for each observation per row, with column names \code{"Fit"} and \code{"Cluster"}.
This \code{data.frame} therefore has \code{nrow(model.data(object)) * nClusters(object)} rows.}
\item{\code{matrix}}{An N-by-K \code{matrix} where each row provides the cluster-specific predictions for the respective observations in \code{newdata}.
Here, \code{N = nrow(newdata)} and \code{K = nClusters(object)}.}
\item{\code{vector}}{A \code{vector} of length \code{nrow(newdata)} with predictions corresponding to the rows of \code{newdata}.}
}

Users can implement support for other prediction formats by defining the \code{transformPredict()} method with other signatures.
}
\section{Example implementation}{

In case we have a custom \code{lcModel} class based on an existing internal model representation with a \code{predict()} function,
we can use \code{transformPredict()} to easily transform the internal model predictions to the right format.
A common output is a \code{matrix} with the cluster-specific predictions.
\preformatted{
predict.lcModelExample <- function(object, newdata) {
  predictionMatrix <- predict(object@model, newdata)
  transformPredict(
    pred = predictionMatrix,
    model = object,
    newdata = newdata
  )
}
}

However, for ease of implementation it is generally advisable to implement \link{predictForCluster} instead of \link{predict.lcModel}.

For a complete and runnable example, see the custom models vignette accessible via \code{vignette("custom", package = "latrend")}.
}

\seealso{
predictForCluster, predict.lcModel
}
