% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/latrend.R
\name{latrend}
\alias{latrend}
\title{Cluster longitudinal data}
\usage{
latrend(
  method,
  data,
  ...,
  envir = NULL,
  verbose = getOption("latrend.verbose")
)
}
\arguments{
\item{method}{The \code{lcMethod} object specifying the longitudinal cluster method to apply.}

\item{data}{The \code{data.frame} or \code{matrix} to which to apply the method.}

\item{...}{Any other arguments to update the \code{lcMethod} definition with.}

\item{envir}{The \code{environment} in which to evaluate the method arguments. Note that this only applies to \code{data} when \code{data} is a \code{call}.}

\item{verbose}{The level of verbosity. Either an object of class \code{Verbose} (see \link[R.utils:Verbose]{R.utils::Verbose} for details),
a \code{logical} indicating whether to show basic computation information,
a \code{numeric} indicating the verbosity level (see \link{Verbose}),
or one of \code{c('info', 'fine', 'finest')}.}
}
\value{
A \code{lcModel} object representing the fitted model.
}
\description{
Cluster longitudinal data
}
\details{
If a seed value is specified in the \code{lcMethod} object or arguments to \code{latrend}, this seed is set using \code{set.seed} prior to the cluster preparation step.
}
\examples{
data(latrendData)
model <- latrend(lcMethodKML("Y", id = "Id", time = "Time"), data = latrendData)

method <- lcMethodKML("Y", id = "Id", time = "Time")
model <- latrend(method, data = latrendData, nClusters = 3)

model <- latrend(method, data = latrendData, nClusters = 3, seed = 1)
}
\seealso{
Other longitudinal cluster fit functions: 
\code{\link{latrendBatch}()},
\code{\link{latrendBoot}()},
\code{\link{latrendCV}()},
\code{\link{latrendRep}()}
}
\concept{longitudinal cluster fit functions}
