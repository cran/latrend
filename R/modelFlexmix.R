#' @include model.R
setClass('lcModelFlexmix', contains = 'lcModel')


#' @export
#' @rdname interface-flexmix
#' @inheritParams fitted.lcModel
fitted.lcModelFlexmix = function(object, ..., clusters = trajectoryAssignments(object)) {
  predNames = paste0('pred_m', 1:nClusters(object))
  predMat = flexmix::fitted(object@model) %>%
    set_colnames(clusterNames(object))
  transformFitted(pred = predMat,
                  model = object,
                  clusters = clusters)
}


#. predictForCluster ####
#' @rdname interface-flexmix
#' @inheritParams predictForCluster
setMethod('predictForCluster', 'lcModelFlexmix', function(
  object, newdata, cluster, what = 'mu', ...)
{
  assert_that(what == 'mu', msg = 'only what="mu" is supported')

  if (nrow(newdata) == 0) {
    return(numeric())
  }

  predOut = flexmix::predict(object@model, newdata = newdata)

  clusIdx = match(cluster, clusterNames(object))

  pred = predOut[[clusIdx]]
  assert_that(ncol(pred) == 1,
    msg = 'unexpected output. the lcModel implementation does not support this model')

  as.numeric(pred)
})


#' @rdname interface-flexmix
setMethod('postprob', 'lcModelFlexmix', function(object, ...) {
  pp = postProbFromObs(object@model@posterior$scaled, make.idRowIndices(object))
  colnames(pp) = clusterNames(object)
  return(pp)
})


#' @export
#' @rdname interface-flexmix
logLik.lcModelFlexmix = function(object, ...) {
  ll = object@model@logLik
  N = nIds(object)
  df = length(coef(object))
  attr(ll, 'nobs') = N
  attr(ll, 'df') = df
  class(ll) = 'logLik'
  ll
}


#' @export
#' @rdname interface-flexmix
coef.lcModelFlexmix = function(object, ...) {
  flexmix::parameters(object@model)
}


#' @rdname interface-flexmix
setMethod('converged', 'lcModelFlexmix', function(object, ...) {
  object@model@converged
})
