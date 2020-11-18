#' @include model.R
setClass('lcModelMixtoolsRM', contains = 'lcModel')


#' @export
#' @importFrom plyr alply
#' @rdname interface-mixtools
#' @inheritParams predict.lcApproxModel
#' @param se Whether to compute the standard error of the prediction.
#' @param ci The confidence interval to compute.
predict.lcModelMixtoolsRM = function(object,
                                     ...,
                                     newdata = NULL,
                                     what = 'mu',
                                     se = TRUE,
                                     ci = c(.025, .975),
                                     approxFun = approx) {
  assert_that(
    is.newdata(newdata),
    what %in% c('mu', 'sigma'),
    is.logical(se),
    is.null(ci) || is.numeric(ci),
    is.function(approxFun)
  )

  if (is.null(newdata)) {
    times = time(object)
    newdata = data.table(Id = rep(ids(object), each = length(times)),
                         Time = times) %>%
      setnames('Id', idVariable(object)) %>%
      setnames('Time', timeVariable(object))
  }
  assert_that(has_name(newdata, timeVariable(object)))

  # compute cluster trajectories
  blocks = object@model$blockid
  assert_that(length(blocks) == length(time(object)))
  respFun = switch(what, mu = mean, sigma = sd)
  if (isTRUE(se)) {
    seFun = function(x)
      sd(x, na.rm = TRUE) / sum(is.finite(x))
  } else {
    seFun = function(x)
      numeric(0)
  }

  comboMat = cbind(
    .Component = seq_len(nClusters(object)) %>% rep(each = length(blocks)),
    .Block = rep(blocks, nClusters(object))
  )

  statMat = plyr::alply(comboMat, 1, function(x) {
    dd = density(
      object@model,
      component = x[1],
      block = x[2],
      scale = FALSE
    )
    c(
      Fit = respFun(dd$x, na.rm = TRUE),
      Se.fit = seFun(dd$x),
      Q = quantile(dd$x, ci)
    )
  }) %>%
    do.call(rbind, .)
  dtStats = cbind(comboMat, statMat) %>% as.data.table()

  newtimes = sort(unique(newdata[[timeVariable(object)]]))
  dtPred = dtStats[, lapply(.SD, function(y)
    approxFun(
      x = time(object),
      y = y,
      xout = newtimes
    )$y), keyby = .Component, .SDcols = -c('.Block')] %>%
    .[, Cluster := factor(.Component,
                          levels = seq_len(nClusters(object)),
                          labels = clusterNames(object))] %>%
    .[, .Component := NULL] %>%
    .[, c(timeVariable(object)) := rep(newtimes, nClusters(object))] %>%
    setcolorder(c('Cluster', timeVariable(object)))

  transformPredict(pred = dtPred,
                   model = object,
                   newdata = newdata)
}


#' @export
#' @rdname interface-mixtools
fitted.lcModelMixtoolsRM = function(object, ..., clusters) {
  predList = predict.lcModelMixtoolsRM(object,
                                       newdata = NULL,
                                       se = FALSE,
                                       ci = NULL)
  transformFitted(predList, model = object, clusters = clusters)
}

#' @rdname interface-mixtools
setMethod('postprob', signature('lcModelMixtoolsRM'), function(object, ...) {
  pp = object@model$posteriors
  colnames(pp) = clusterNames(object)
  return(pp)
})


#' @export
#' @rdname interface-mixtools
logLik.lcModelMixtoolsRM = function(object, ...) {
  ll = object@model$loglik
  attr(ll, 'nobs') = nIds(object)
  attr(ll, 'df') = length(coef(object)) + 1
  class(ll) = 'logLik'
  return(ll)
}

#' @rdname interface-mixtools
setMethod('converged', signature('lcModelMixtoolsRM'), function(object, ...) {
  TRUE
})