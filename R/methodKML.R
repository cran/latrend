#' @include method.R

#' @name interface-kml
#' @rdname interface-kml
#' @title kml interface
#' @seealso [lcMethodKML] \link[kml]{kml}
#' @keywords internal
NULL

setClass('lcMethodKML', contains = 'lcMatrixMethod')

#' @export
#' @title Specify a longitudinal k-means (KML) method
#' @param response The name of the response variable.
#' @param time The name of the time variable.
#' @param id The name of the trajectory identifier variable.
#' @param nClusters The number of clusters to estimate.
#' @param ... Arguments passed to [kml::parALGO] and [kml::kml].
#' The following external arguments are ignored: object, nbClusters, parAlgo, toPlot, saveFreq
#' @examples
#' library(kml)
#' data(latrendData)
#' method <- lcMethodKML("Y", id = "Id", time = "Time", nClusters = 3)
#' model <- latrend(method, latrendData)
#' @references
#' \insertRef{genolini2015kml}{latrend}
#' @family lcMethod implementations
lcMethodKML = function(response,
                       time = getOption('latrend.time'),
                       id = getOption('latrend.id'),
                       nClusters = 2,
                       ...) {
  lcMethod.call(
    'lcMethodKML',
    call = match.call.defaults(),
    defaults = c(kml::kml, kml::parALGO),
    excludeArgs = c('object', 'nbClusters', 'parAlgo', 'toPlot', 'saveFreq')
  )
}

#' @rdname interface-kml
setMethod('getName', signature('lcMethodKML'), function(object) 'longitudinal k-means (KML)')

#' @rdname interface-kml
setMethod('getShortName', signature('lcMethodKML'), function(object) 'kml')

#' @rdname interface-kml
#' @inheritParams preFit
setMethod('preFit', signature('lcMethodKML'), function(method, data, envir, verbose, ...) {
  e = callNextMethod()

  valueColumn = responseVariable(method)

  # Model specification
  cat(verbose, 'Creating clusterLongData object...', level = verboseLevels$finest)

  parRefArgs = list(saveFreq = Inf, scale = FALSE)
  parArgs = modifyList(parRefArgs, as.list(method, args = kml::parALGO), keep.null = TRUE)
  e$par = do.call(kml::parALGO, parArgs)

  e$cld = kml::clusterLongData(
    traj = e$dataMat,
    idAll = rownames(e$dataMat),
    time = sort(unique(data[[timeVariable(method)]]))
  )
  return(e)
})

#' @rdname interface-kml
#' @inheritParams fit
setMethod('fit', signature('lcMethodKML'), function(method, data, envir, verbose, ...) {
  cld = envir$cld

  # Helper variables
  valueColumn = responseVariable(method)
  suppressFun = ifelse(as.logical(verbose), force, capture.output)

  cat(verbose, 'Running kml()...', level = verboseLevels$finest)
  suppressFun(
    # note that slowKML throws an error for nbClusters=1
    kml::kml(
      cld,
      nbClusters = method$nClusters,
      nbRedrawing = method$nbRedrawing,
      toPlot = 'none',
      parAlgo = envir$par
    )
  )

  new(
    'lcModelKML',
    method = method,
    data = data,
    model = cld,
    clusterNames = make.clusterNames(method$nClusters)
  )
})
