#' @include model.R
setClass('lcModelLMKM',
         representation(coefNames = 'character'),
         contains = 'lcModel')

#' @export
#' @rdname interface-featureBased
#' @param cluster The cluster name.
coef.lcModelLMKM = function(object, ..., cluster = NULL) {
  coefmat = t(object@model$centers)
  colnames(coefmat) = clusterNames(object)
  rownames(coefmat) = object@coefNames

  if (is.null(cluster)) {
    return(coefmat)
  } else {
    assert_that(is.count(cluster) ||
                  is.character(cluster) && cluster %in% clusterNames(object))
    coefmat[, cluster, drop = FALSE]
  }
}

#. converged ####
#' @rdname interface-featureBased
setMethod('converged', signature('lcModelLMKM'), function(object, ...) {
  if (nClusters(object) == 1) {
    TRUE
  }
  else {
    !object@model$ifault
  }
})

#. postprob ####
#' @rdname interface-featureBased
setMethod('postprob', signature('lcModelLMKM'), function(object, ...) {
  k = nrow(object@model$centers)
  postprobFromAssignments(object@model$cluster, k)
})

#' @export
#' @rdname interface-featureBased
#' @inheritParams predict.lcModel
#' @inheritDotParams stats::predict.lm
predict.lcModelLMKM = function(object, ...,
                               newdata = NULL,
                               what = 'mu') {
  assert_that(is.newdata(newdata))
  assert_that(what == 'mu')

  if (is.null(newdata)) {
    newdata = model.data(object) %>%
      subset(select = getCovariates(formula(object)))
  }
  newdata = as.data.table(newdata)

  if (nrow(newdata) == 0) {
    # predict.lm cannot handle empty data.frame, so return early
    return(transformPredict(
      pred = NULL,
      model = object,
      newdata = newdata
    ))
  }

  # create ref lm
  method = getLcMethod(object)
  id = idVariable(method)
  lmArgs = as.list(method, args = lm)
  data = model.data(object)
  refdata = data[get(id) == first(get(id))][1:min(10, .N)]
  refmod = do.call(lm, c(lmArgs, data = list(refdata)))

  # construct lm per cluster
  clusmods = lapply(clusterNames(object), function(clusName) {
    clusmod = refmod
    clusmod$coefficients = coef(object, cluster = clusName)
    return(clusmod)
  })

  predfun = function(clusmod, clusdata) {
    out = predict(clusmod, newdata = clusdata, ...)
    if (is.numeric(out)) {
      list(fit = out)
    } else {
      out
    }
  }

  if (has_name(newdata, 'Cluster')) {
    clusdataList = split(newdata, by = 'Cluster', sorted = TRUE)
  } else {
    clusdataList = list(newdata)
  }

  dtpred = mapply(predfun, clusmods, clusdataList, SIMPLIFY = FALSE) %>%
    rbindlist(idcol = 'Cluster') %>%
    .[, Cluster := factor(Cluster, labels = clusterNames(object))]
  setnames(dtpred, 'fit', 'Fit')
  transformPredict(pred = dtpred,
                   model = object,
                   newdata = newdata)
}