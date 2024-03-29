% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/model.R, R/metricsInternal.R,
%   R/models.R
\name{metric}
\alias{metric}
\alias{metric,lcModel-method}
\alias{internalMetric}
\alias{metric,list-method}
\alias{metric,lcModels-method}
\title{Compute internal model metric(s)}
\usage{
metric(object, name = getOption("latrend.metric", c("WRSS", "APPA.mean")), ...)

\S4method{metric}{lcModel}(object, name = getOption("latrend.metric", c("WRSS", "APPA.mean")), ...)

\S4method{metric}{list}(object, name, drop = TRUE)

\S4method{metric}{lcModels}(object, name, drop = TRUE)
}
\arguments{
\item{object}{The \code{lcModel}, \code{lcModels}, or \code{list} of \code{lcModel} objects to compute the metrics for.}

\item{name}{The name(s) of the metric(s) to compute. If no names are given, the names specified in the \code{latrend.metric} option (WRSS, APPA, AIC, BIC) are used.}

\item{...}{Additional arguments.}

\item{drop}{Whether to return a \verb{numeric vector} instead of a \code{data.frame}
in case of a single metric.}
}
\value{
For \code{metric(lcModel)}: A named \code{numeric} vector with the computed model metrics.

For \code{metric(list)}: A \code{data.frame} with a metric per column.

For \code{metric(lcModels)}: A \code{data.frame} with a metric per column.
}
\description{
Compute one or more internal metrics for the given \code{lcModel} object.

Note that there are many metrics available, and there exists no metric that works best in all scenarios.
It is recommended to carefully consider which metric is most appropriate for your use case.

Recommended overview papers:
\itemize{
\item \insertCite{arbelaitz2013extensive;textual}{latrend} provide an extensive overview validity indices for cluster algorithms.
\item \insertCite{vandernest2020overview;textual}{latrend} provide an overview of metrics for mixture models (GBTM, GMM); primarily likelihood-based or posterior probability-based metrics.
\item \insertCite{henson2007detecting;textual}{latrend} provide an overview of likelihood-based metrics for mixture models.
}

Call \code{\link[=getInternalMetricNames]{getInternalMetricNames()}} to retrieve the names of the defined internal metrics.

See the \emph{Details} section below for a list of supported metrics.
}
\section{Supported internal metrics}{
\tabular{lll}{
   \strong{Metric name} \tab \strong{Description} \tab \strong{Function / Reference} \cr
   \code{AIC} \tab \href{https://en.wikipedia.org/wiki/Akaike_information_criterion}{Akaike information criterion}. A goodness-of-fit estimator that adjusts for model complexity (i.e., the number of parameters). Only available for models that support the computation of the model log-likelihood through \link{logLik}. \tab \code{\link[stats:AIC]{stats::AIC()}}, \insertCite{akaike1974new}{latrend} \cr
   \code{APPA.mean} \tab Mean of the average posterior probability of assignment (APPA) across clusters. A measure of the precision of the trajectory classifications. A score of 1 indicates perfect classification. \tab \code{\link[=APPA]{APPA()}}, \insertCite{nagin2005group}{latrend} \cr
   \code{APPA.min} \tab Lowest APPA among the clusters \tab \code{\link[=APPA]{APPA()}}, \insertCite{nagin2005group}{latrend} \cr
   \code{ASW} \tab Average \href{https://en.wikipedia.org/wiki/Silhouette_(clustering)}{silhouette} width based on the Euclidean distance \tab \insertCite{rousseeuw1987silhouettes}{latrend} \cr
   \code{BIC} \tab \href{https://en.wikipedia.org/wiki/Bayesian_information_criterion}{Bayesian information criterion}. A goodness-of-fit estimator that corrects for the degrees of freedom (i.e., the number of parameters) and sample size. Only available for models that support the computation of the model log-likelihood through \link{logLik}. \tab \code{\link[stats:AIC]{stats::BIC()}}, \insertCite{schwarz1978estimating}{latrend} \cr
   \code{CAIC} \tab Consistent Akaike information criterion \tab \insertCite{bozdogan1987model}{latrend} \cr
   \code{CLC} \tab Classification likelihood criterion \tab \insertCite{mclachlan2000finite}{latrend} \cr
   \code{converged} \tab Whether the model converged during estimation \tab \code{\link[=converged]{converged()}} \cr
   \code{deviance} \tab The model \href{https://en.wikipedia.org/wiki/Deviance_(statistics)}{deviance} \tab \code{\link[stats:deviance]{stats::deviance()}} \cr
   \code{Dunn} \tab The \href{https://en.wikipedia.org/wiki/Dunn_index}{Dunn index} \tab \insertCite{dunn1974well}{latrend} \cr
   \code{entropy} \tab Entropy of the posterior probabilities \tab  \cr
   \code{estimationTime} \tab The time needed for fitting the model \tab \code{\link[=estimationTime]{estimationTime()}} \cr
   \code{ED} \tab \href{https://en.wikipedia.org/wiki/Euclidean_distance}{Euclidean distance} between the cluster trajectories and the assigned observed trajectories \tab  \cr
   \code{ED.fit} \tab Euclidean distance between the cluster trajectories and the assigned fitted trajectories \tab  \cr
   \code{ICL.BIC} \tab Integrated classification likelihood (ICL) approximated using the BIC \tab \insertCite{biernacki2000assessing}{latrend} \cr
   \code{logLik} \tab Model log-\href{https://en.wikipedia.org/wiki/Likelihood_function}{likelihood} \tab \code{\link[stats:logLik]{stats::logLik()}} \cr
   \code{MAE} \tab \href{https://en.wikipedia.org/wiki/Mean_absolute_error}{Mean absolute error} of the fitted trajectories (assigned to the most likely respective cluster) to the observed trajectories \tab  \cr
   \code{Mahalanobis} \tab \href{https://en.wikipedia.org/wiki/Mahalanobis_distance}{Mahalanobis distance} between the cluster trajectories and the assigned observed trajectories \tab \insertCite{mahalanobis1936generalized}{latrend} \cr
   \code{MSE} \tab \href{https://en.wikipedia.org/wiki/Mean_squared_error}{Mean squared error} of the fitted trajectories (assigned to the most likely respective cluster) to the observed trajectories \tab  \cr
   \code{relativeEntropy}, \code{RE} \tab A measure of the precision of the trajectory classification. A value of 1 indicates perfect classification, whereas a value of 0 indicates a non-informative uniform classification. It is the normalized version of \code{entropy}, scaled between [0, 1]. \tab \insertCite{ramaswamy1993empirical}{latrend}, \insertCite{muthen2004latent}{latrend} \cr
   \code{RMSE} \tab \href{https://en.wikipedia.org/wiki/Root-mean-square_deviation}{Root mean squared error} of the fitted trajectories (assigned to the most likely respective cluster) to the observed trajectories \tab  \cr
   \code{RSS} \tab \href{https://en.wikipedia.org/wiki/Residual_sum_of_squares}{Residual sum of squares} under most likely cluster allocation \tab  \cr
   \code{scaledEntropy} \tab See \code{relativeEntropy} \tab  \cr
   \code{sigma} \tab The residual standard deviation \tab \code{\link[stats:sigma]{stats::sigma()}} \cr
   \code{ssBIC} \tab Sample-size adjusted BIC \tab \insertCite{sclove1987application}{latrend} \cr
   \code{SED} \tab Standardized Euclidean distance between the cluster trajectories and the assigned observed trajectories \tab  \cr
   \code{SED.fit} \tab The cluster-weighted standardized Euclidean distance between the cluster trajectories and the assigned fitted trajectories \tab  \cr
   \code{WMAE} \tab \code{MAE} weighted by cluster-assignment probability \tab  \cr
   \code{WMSE} \tab \code{MSE} weighted by cluster-assignment probability \tab  \cr
   \code{WRMSE} \tab \code{RMSE} weighted by cluster-assignment probability \tab  \cr
   \code{WRSS} \tab \code{RSS} weighted by cluster-assignment probability \tab  \cr
}
}

\section{Implementation}{

See the documentation of the \code{\link[=defineInternalMetric]{defineInternalMetric()}} function for details on how to define your own metrics.
}

\examples{
data(latrendData)
method <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time")
model <- latrend(method, latrendData)
metric(model, "WMAE")

if (require("clusterCrit")) {
  metric(model, c("WMAE", "Dunn"))
}
}
\references{
\insertAllCited{}
}
\seealso{
\link{externalMetric} \link{min.lcModels} \link{max.lcModels}

Other metric functions: 
\code{\link{defineExternalMetric}()},
\code{\link{defineInternalMetric}()},
\code{\link{externalMetric}()},
\code{\link{getExternalMetricDefinition}()},
\code{\link{getExternalMetricNames}()},
\code{\link{getInternalMetricDefinition}()},
\code{\link{getInternalMetricNames}()}

Other lcModel functions: 
\code{\link{clusterNames}()},
\code{\link{clusterProportions}()},
\code{\link{clusterSizes}()},
\code{\link{clusterTrajectories}()},
\code{\link{coef.lcModel}()},
\code{\link{converged}()},
\code{\link{deviance.lcModel}()},
\code{\link{df.residual.lcModel}()},
\code{\link{estimationTime}()},
\code{\link{externalMetric}()},
\code{\link{fitted.lcModel}()},
\code{\link{fittedTrajectories}()},
\code{\link{getCall.lcModel}()},
\code{\link{getLcMethod}()},
\code{\link{ids}()},
\code{\link{lcModel-class}},
\code{\link{model.frame.lcModel}()},
\code{\link{nClusters}()},
\code{\link{nIds}()},
\code{\link{nobs.lcModel}()},
\code{\link{plot-lcModel-method}},
\code{\link{plotClusterTrajectories}()},
\code{\link{plotFittedTrajectories}()},
\code{\link{postprob}()},
\code{\link{predict.lcModel}()},
\code{\link{predictAssignments}()},
\code{\link{predictForCluster}()},
\code{\link{predictPostprob}()},
\code{\link{qqPlot}()},
\code{\link{residuals.lcModel}()},
\code{\link{sigma.lcModel}()},
\code{\link{strip}()},
\code{\link{time.lcModel}()},
\code{\link{trajectoryAssignments}()}
}
\concept{lcModel functions}
\concept{metric functions}
