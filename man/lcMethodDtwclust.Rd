% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodDtwclust.R
\name{lcMethodDtwclust}
\alias{lcMethodDtwclust}
\title{Specify time series clustering via dtwclust}
\usage{
lcMethodDtwclust(
  response,
  time = getOption("latrend.time"),
  id = getOption("latrend.id"),
  nClusters = 2,
  ...
)
}
\arguments{
\item{response}{The name of the response variable.}

\item{time}{The name of the time variable.}

\item{id}{The name of the trajectory identifier variable.}

\item{nClusters}{Number of clusters.}

\item{...}{Arguments passed to \link[dtwclust:tsclust]{dtwclust::tsclust}.
The following arguments are ignored: series, k, trace.}
}
\description{
Specify time series clustering via dtwclust
}
\examples{
data(latrendData)

if (require("dtwclust")) {
  method <- lcMethodDtwclust("Y", id = "Id", time = "Time", nClusters = 3)
  # reduced sample size to lower runtime
  model <- latrend(method, latrendData[1:500, ])
}
}
\references{
\insertRef{sardaespinosa2019time}{latrend}
}
\seealso{
Other lcMethod implementations: 
\code{\link{getArgumentDefaults}()},
\code{\link{getArgumentExclusions}()},
\code{\link{lcMethod-class}},
\code{\link{lcMethodAkmedoids}},
\code{\link{lcMethodCrimCV}},
\code{\link{lcMethodFeature}},
\code{\link{lcMethodFunFEM}},
\code{\link{lcMethodFunction}},
\code{\link{lcMethodGCKM}},
\code{\link{lcMethodKML}},
\code{\link{lcMethodLMKM}},
\code{\link{lcMethodLcmmGBTM}},
\code{\link{lcMethodLcmmGMM}},
\code{\link{lcMethodMclustLLPA}},
\code{\link{lcMethodMixAK_GLMM}},
\code{\link{lcMethodMixtoolsGMM}},
\code{\link{lcMethodMixtoolsNPRM}},
\code{\link{lcMethodRandom}},
\code{\link{lcMethodStratify}}
}
\concept{lcMethod implementations}
