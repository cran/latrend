## ----setup, include = FALSE---------------------------------------------------
library(latrend)
library(ggplot2)
set.seed(1)
knitr::opts_chunk$set(
  cache = TRUE,
  collapse = TRUE,
  fig.width = 7,
  fig.align = "center",
  fig.topcaption = TRUE,
  comment = "#>"
)

## ---- results='hide',message=FALSE,warning=FALSE------------------------------
library(latrend)
data(latrendData)

## -----------------------------------------------------------------------------
head(latrendData)

## -----------------------------------------------------------------------------
options(latrend.id = "Id", latrend.time = "Time")

## ---- fig.asp = .6, fig.cap='Visualizing the trajectories of the `latrend` dataset.'----
plotTrajectories(latrendData, response = "Y")

## -----------------------------------------------------------------------------
kmlMethod <- lcMethodKML(response = "Y", nClusters = 2, nbRedrawing = 1)

kmlMethod

## -----------------------------------------------------------------------------
kmlModel <- latrend(kmlMethod, data = latrendData)

## -----------------------------------------------------------------------------
kmlModel

## -----------------------------------------------------------------------------
kmlMethods <- lcMethods(kmlMethod, nClusters = 1:8)

as.data.frame(kmlMethods)

## -----------------------------------------------------------------------------
kmlModels <- latrendBatch(kmlMethods, data = latrendData, verbose = FALSE)

kmlModels

## ----warning=FALSE, fig.cap = 'Elbow plots of three relevant cluster metrics across the fitted models.'----
plotMetric(kmlModels, c("logLik", "BIC", "WMAE"))

## -----------------------------------------------------------------------------
kmlModel4 <- subset(kmlModels, nClusters == 4, drop = TRUE)

kmlModel4

## ---- fig.cap = 'Cluster trajectories for KML model with 4 clusters.'---------
plotClusterTrajectories(kmlModel4)

## ---- fig.asp = .8, fig.cap = 'Cluster trajectories for KML model with 4 clusters, along with the assigned trajectories.'----
plot(kmlModel4)

## -----------------------------------------------------------------------------
getInternalMetricNames()

## -----------------------------------------------------------------------------
metric(kmlModel, c("APPA", "WRSS", "WMAE"))

## ---- fig.width = 5, fig.cap = 'QQ-plot of the selected KML model.'-----------
qqPlot(kmlModel4)

## ---- fig.asp = .8, fig.cap = 'Cluster-specific detrended QQ-plot for the selected KML model.'----
qqPlot(kmlModel4, byCluster = TRUE, detrend = TRUE)

## -----------------------------------------------------------------------------
library(splines)
gbtmMethod <- lcMethodLcmmGBTM(Y ~ bs(Time) * CLUSTER + (1 | Id))

gbtmMethod

## -----------------------------------------------------------------------------
gbtmMethods <- lcMethods(gbtmMethod, nClusters = 1:5)

gbtmModels <- latrendBatch(gbtmMethods, data = latrendData, verbose = FALSE)

## ----warning=FALSE, fig.cap = 'Three cluster metrics for each of the GBTMs.'----
plotMetric(gbtmModels, c("logLik", "BIC", "WMAE"))

## -----------------------------------------------------------------------------
bestGbtmModel <- subset(gbtmModels, nClusters == 3, drop=TRUE)
plot(bestGbtmModel)

## -----------------------------------------------------------------------------
gmmMethod <- lcMethodLcmmGMM(Y ~ poly(Time, 2, raw = TRUE) * CLUSTER + (1 | Id), idiag = TRUE)

gmmMethod

## -----------------------------------------------------------------------------
gmmMethods <- lcMethods(gmmMethod, nClusters = 1:5)

gmmModels <- latrendBatch(gmmMethods, latrendData, verbose = FALSE)

## ----warning=FALSE, fig.cap = 'Three cluster metrics for each of the GMMs.'----
plotMetric(gmmModels, c("logLik", "BIC", "WMAE"))

## -----------------------------------------------------------------------------
bestGmmModel <- subset(gmmModels, nClusters == 3, drop=TRUE)

## ---- fig.cap = 'Cluster trajectories of the selected GMM, including the assigned trajectories.'----
plot(bestGmmModel)

## ---- fig.width = 5, fig.cap = 'Detrended QQ plot for the best GMM.'----------
qqPlot(bestGmmModel, detrend = TRUE)

## -----------------------------------------------------------------------------
BIC(bestGbtmModel, bestGmmModel)

## -----------------------------------------------------------------------------
metric(list(bestGbtmModel, bestGmmModel), 'WMAE')

## -----------------------------------------------------------------------------
externalMetric(bestGbtmModel, bestGmmModel, 'WMMAE')

## -----------------------------------------------------------------------------
externalMetric(bestGbtmModel, bestGmmModel, 'adjustedRand')

## ---- fig.width = 5, fig.cap = 'Non-parametric estimates of the cluster trajectories based on the reference assignments.'----
plotClusterTrajectories(latrendData, response = "Y", cluster = "Class")

## -----------------------------------------------------------------------------
refTrajAssigns <- aggregate(Class ~ Id, data = latrendData, FUN = data.table::first)
refModel <- lcModelPartition(data = latrendData, response = "Y", trajectoryAssignments = refTrajAssigns$Class)
refModel

## ---- fig.cap = 'Cluster trajectories of the reference model.'----------------
plot(refModel)

## -----------------------------------------------------------------------------
getExternalMetricNames()

## -----------------------------------------------------------------------------
externalMetric(bestGmmModel, refModel, "adjustedRand")

