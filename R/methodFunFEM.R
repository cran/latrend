#' @include method.R

#' @name interface-funFEM
#' @rdname interface-funFEM
#' @title funFEM interface
#' @seealso [lcMethodFunFEM] \link[funFEM]{funFEM-package}
#' @keywords internal
NULL

setClass('lcMethodFunFEM', contains = 'lcMatrixMethod')

#' @export
#' @title Specify a FunFEM method
#' @inheritParams lcMatrixMethod-class
#' @inheritParams lcMethodKML
#' @param basis The basis function. By default, a 3rd-order B-spline with 10 breaks is used.
#' @param ... Arguments passed to [funFEM::funFEM].
#' The following external arguments are ignored: fd, K, disp, graph.
#' @examples
#' library(funFEM)
#' library(fda)
#' data(latrendData)
#' method <- lcMethodFunFEM("Y", id = "Id", time = "Time", nClusters = 3)
#' model <- latrend(method, latrendData)
#'
#' method <- lcMethodFunFEM("Y",
#'    basis = function(time) {
#'       create.bspline.basis(time,
#'         nbasis = 10, norder = 4)
#' })
#' @family lcMethod implementations
#' @references
#' \insertRef{bouveyron2015funfem}{latrend}
lcMethodFunFEM = function(response,
                          time = getOption('latrend.time'),
                          id = getOption('latrend.id'),
                          nClusters = 2,
                          basis = function(time) fda::create.bspline.basis(time, nbasis = 10, norder = 4),
                          ...) {
  lcMethod.call(
    'lcMethodFunFEM',
    call = match.call.defaults(),
    defaults = funFEM::funFEM,
    excludeArgs = c('fd', 'K', 'disp', 'graph')
  )
}

#' @rdname interface-funFEM
setMethod('getName', signature('lcMethodFunFEM'), function(object) 'functional subspace clustering with FunFEM')

#' @rdname interface-funFEM
setMethod('getShortName', signature('lcMethodFunFEM'), function(object) 'funfem')

#' @rdname interface-funFEM
setMethod('preFit', signature('lcMethodFunFEM'), function(method, data, envir, verbose, ...) {
  requireNamespace('fda')
  requireNamespace('funFEM')

  e = callNextMethod()
  e$basis = method$basis(range(e$times))
  e$fd = fda::smooth.basis(e$times, t(e$dataMat), e$basis)$fd

  return(e)
})

#' @rdname interface-funFEM
#' @inheritParams fit
setMethod('fit', signature('lcMethodFunFEM'), function(method, data, envir, verbose, ...) {
  args = as.list(method, args = funFEM::funFEM)
  args$fd = envir$fd
  args$K = method$nClusters
  args$disp = FALSE
  args$graph = FALSE

  # Helper variables
  suppressFun = ifelse(as.logical(verbose), force, capture.output)

  suppressFun({
    model = do.call(funFEM::funFEM, args)
  })
  model$basis = envir$basis
  model$fd = envir$fd

  new(
    'lcModelFunFEM',
    method = method,
    data = data,
    model = model,
    clusterNames = make.clusterNames(method$nClusters)
  )
})
