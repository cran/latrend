% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method.R
\name{lcMethod-estimation}
\alias{lcMethod-estimation}
\alias{latrend-procedure}
\alias{lcMethod-steps}
\title{Longitudinal cluster method (\code{lcMethod}) estimation procedure}
\description{
Each longitudinal cluster method represented by a \link[=lcMethod-class]{lcMethod class} implements a series of standardized steps that produce the estimated method as its output.
These steps, as part of the estimation procedure, are executed by the \code{\link[=latrend]{latrend()}} function and other functions prefixed by \emph{"latrend"} (e.g., \code{\link[=latrendRep]{latrendRep()}}, \code{\link[=latrendBoot]{latrendBoot()}}, \code{\link[=latrendCV]{latrendCV()}}).
}
\section{Estimation procedure}{

The steps for estimating a \code{lcMethod} object are defined and executed as follows:
\enumerate{
\item \code{\link[=compose]{compose()}}: Evaluate and finalize the method argument values.
\item \code{\link[=validate]{validate()}}: Check the validity of the method argument values in relation to the dataset.
\item \code{\link[=prepareData]{prepareData()}}: Process the training data for fitting.
\item \code{\link[=preFit]{preFit()}}: Prepare environment for estimation, independent of training data.
\item \code{\link[=fit]{fit()}}: Estimate the specified method on the training data, outputting an object inheriting from \code{lcModel}.
\item \code{\link[=postFit]{postFit()}}: Post-process the outputted \code{lcModel} object.
}

The result of the fitting procedure is an \linkS4class{lcModel} object that inherits from the \code{lcModel} class.
}

\examples{
data(latrendData)
method <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time")
model <- latrend(method, data = latrendData)
summary(model)
}
\seealso{
\linkS4class{lcMethod} \link{latrend}
}
