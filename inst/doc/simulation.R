## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = TRUE
)

## ----setup, include = FALSE---------------------------------------------------
library(ggplot2)
library(magrittr)
library(latrend)
library(simTool)
knitr::opts_chunk$set(
  cache = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(latrend)
options(latrend.verbose = FALSE)

## -----------------------------------------------------------------------------
dataGen <- function(numTraj, ..., data.seed) {
  latrend::generateLongData(
    sizes = c(floor(numTraj * .4), ceiling(numTraj * .6)),
    fixed = Y ~ 1 + Time,
    cluster = ~ 1 + Time,
    random = ~ 1,
    id = "Traj",
    clusterCoefs = cbind(c(0, 0), c(1, -1)),
    seed = data.seed,
    ...
  )
}

## -----------------------------------------------------------------------------
options(latrend.id = "Traj", latrend.time = "Time")

## -----------------------------------------------------------------------------
exampleData <- dataGen(numTraj = 200, randomScale = .1, data.seed = 1)
plotTrajectories(exampleData, response = "Y")

## -----------------------------------------------------------------------------
dataGrid <- simTool::expand_tibble(
  fun = "dataGen",
  numTraj = c(50, 250),
  randomScales = c(.1, .5),
  data.seed = 1:2
)

head(dataGrid)

## -----------------------------------------------------------------------------
kmlMethodGrid <- simTool::expand_tibble(
  proc = "latrend",
  method = "lcMethodKML",
  nClusters = 1:3,
  seed = 1,
  response = "Y"
)

head(kmlMethodGrid)

## -----------------------------------------------------------------------------
fitGCKM <- function(type, ...) {
  form <- switch(type,
    constant = Y ~ 1 + Time + (1 | Traj),
    linear = Y ~ 1 + Time + (Time | Traj)
  )
  
  latrend(..., formula = form)
}

## -----------------------------------------------------------------------------
gckmMethodGrid <- simTool::expand_tibble(
  proc = "fitGCKM",
  method = "lcMethodGCKM",
  type = c("constant", "linear"),
  nClusters = 1:3,
  seed = 1:2
)

## -----------------------------------------------------------------------------
methodGrid <- dplyr::bind_rows(kmlMethodGrid, gckmMethodGrid)
head(methodGrid)

## -----------------------------------------------------------------------------
analyzeModel <- function(model) {
  data <- model.data(model)
  refModel <- lcModelPartition(data, response = "Y", trajectoryAssignments = "Class")
  
  tibble::tibble(
    BIC = BIC(model),
    APPA = metric(model, "APPA"),
    WMAE = metric(model, "WMAE"),
    ARI = externalMetric(model, refModel, "adjustedRand")
  )
}

## -----------------------------------------------------------------------------
result <- simTool::eval_tibbles(
  data_grid = dataGrid, 
  proc_grid = methodGrid,
  post_analyze = analyzeModel
)

## -----------------------------------------------------------------------------
result

## -----------------------------------------------------------------------------
library(data.table)
resultsTable <- as.data.table(result$simulation)

## -----------------------------------------------------------------------------
resultsTable[, .(K = nClusters[which.min(BIC)]), keyby = .(numTraj, randomScales, data.seed, method)]

## -----------------------------------------------------------------------------
resultsTable[nClusters > 1, .(ARI = mean(ARI)), keyby = .(nClusters, numTraj, randomScales, method)]

## -----------------------------------------------------------------------------
resultsTable[nClusters > 1, .(ARI = mean(ARI)), keyby = .(randomScales, nClusters, method)]

## -----------------------------------------------------------------------------
resultsTable[randomScales == .1, .(WMAE = mean(WMAE)), keyby = .(nClusters, method)]

## -----------------------------------------------------------------------------
resultsTable[, .(WMAE = mean(WMAE)), keyby = .(randomScales, nClusters)]

