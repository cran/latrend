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
options(latrend.verbose = FALSE)

## ---- results='hide',message=FALSE,warning=FALSE------------------------------
library(latrend)

## -----------------------------------------------------------------------------
data(latrendData)
head(latrendData)

## -----------------------------------------------------------------------------
options(latrend.id = "Id", latrend.time = "Time")

## -----------------------------------------------------------------------------
kml <- lcMethodKML("Y", nClusters = 3, nbRedrawing = 5)
kml

## -----------------------------------------------------------------------------
repModels <- latrendRep(kml, data = latrendData, .rep=10)
print(repModels, excludeShared = FALSE)

## -----------------------------------------------------------------------------
repSelfMetrics <- metric(repModels, name = c("BIC", "WMAE", "APPA"))
head(repSelfMetrics)

## -----------------------------------------------------------------------------
summary(repSelfMetrics[, "WMAE"])

## -----------------------------------------------------------------------------
bestRepModel <- min(repModels, "BIC")
externalMetric(repModels, bestRepModel, name = "adjustedRand")

## -----------------------------------------------------------------------------
simMat <- externalMetric(repModels, name = "adjustedRand")
round(simMat, 2)

## -----------------------------------------------------------------------------
summary(simMat)

## -----------------------------------------------------------------------------
bootModels <- latrendBoot(kml, data = latrendData, samples = 10)
bootModels

## -----------------------------------------------------------------------------
bootMetrics <- metric(bootModels, name = c("BIC", "WMAE", "APPA"))
bootMetrics

## -----------------------------------------------------------------------------
colMeans(bootMetrics)

## -----------------------------------------------------------------------------
apply(bootMetrics, 2, sd)

## -----------------------------------------------------------------------------
trainModels <- latrendCV(kml, data = latrendData, folds = 10, seed = 1)
trainModels

## -----------------------------------------------------------------------------
dataFolds <- createTrainDataFolds(latrendData, folds = 10)
foldModels <- latrendBatch(kml, data = dataFolds)
foldModels

## -----------------------------------------------------------------------------
testDataFolds <- createTestDataFolds(latrendData, dataFolds)

