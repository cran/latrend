% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{latrend-methods}
\alias{latrend-methods}
\title{Supported methods for longitudinal clustering}
\description{
This page provides an overview of the currently supported methods for longitudinal clustering.
For general recommendations on which method to apply to your dataset, \link[=latrend-approaches]{see here}.
}
\section{Supported methods}{
\tabular{lll}{
   \strong{Method} \tab \strong{Description} \tab \strong{Source} \cr
   \link{lcMethodAkmedoids} \tab Anchored \emph{k}-medoids \insertCite{adepeju2020akmedoids}{latrend} \tab \code{akmedoids} \cr
   \link{lcMethodCrimCV} \tab Group-based trajectory modeling of count data \insertCite{nielsen2018crimcv}{latrend} \tab \code{crimCV} \cr
   \link{lcMethodDtwclust} \tab Methods for distance-based clustering, including dynamic time warping \insertCite{sardaespinosa2019time}{latrend} \tab \code{dtwclust} \cr
   \link{lcMethodFeature} \tab Feature-based clustering \tab  \cr
   \link{lcMethodFlexmix} \tab Interface to the FlexMix framework \insertCite{gruen2008flexmix}{latrend} \tab \code{flexmix} \cr
   \link{lcMethodFlexmixGBTM} \tab Group-based trajectory modeling \tab \code{flexmix} \cr
   \link{lcMethodFunFEM} \tab Model-based clustering using funFEM \insertCite{bouveyron2015funfem}{latrend} \tab \code{funFEM} \cr
   \link{lcMethodGCKM} \tab Growth-curve modeling and \emph{k}-means \tab \code{lme4} \cr
   \link{lcMethodKML} \tab Longitudinal \emph{k}-means \insertCite{genolini2015kml}{latrend} \tab \code{kml} \cr
   \link{lcMethodLcmmGBTM} \tab Group-based trajectory modeling \insertCite{proustlima2017estimation}{latrend} \tab \code{lcmm} \cr
   \link{lcMethodLcmmGMM} \tab Growth mixture modeling \insertCite{proustlima2017estimation}{latrend} \tab \code{lcmm} \cr
   \link{lcMethodLMKM} \tab Feature-based clustering using linear regression and \emph{k}-means \tab  \cr
   \link{lcMethodMclustLLPA} \tab Longitudinal latent profile analysis \insertCite{scrucca2016mclust}{latrend} \tab \code{mclust} \cr
   \link{lcMethodMixAK_GLMM} \tab Mixture of generalized linear mixed models \tab \code{mixAK} \cr
   \link{lcMethodMixtoolsGMM} \tab Growth mixture modeling \tab \code{mixtools} \cr
   \link{lcMethodMixtoolsNPRM} \tab Non-parametric repeated measures clustering \insertCite{benaglia2009mixtools}{latrend} \tab \code{mixtools} \cr
   \link{lcMethodMixTVEM} \tab Mixture of time-varying effects models \tab  \cr
   \link{lcMethodRandom} \tab Random partitioning \tab  \cr
   \link{lcMethodStratify} \tab Stratification rule \tab  \cr
}


In addition, the functionality of any method can be extended via \link[=lcMetaMethods]{meta methods}.
This is used for extending the estimation procedure of a method, such as \link[=lcFitRep]{repeated fitting} and selecting the best result, or \link[=lcFitConverged]{fitting until convergence}.

It is strongly encouraged to \link[=latrend-metrics]{evaluate and compare} several candidate methods in order to identify the most suitable method.
}

\examples{
data(latrendData)
method <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time")
model <- latrend(method, data = latrendData)
}
\references{
\insertAllCited{}
}
\seealso{
\link{latrend-approaches} \link{latrend-estimation} \link{latrend-metrics}
}
