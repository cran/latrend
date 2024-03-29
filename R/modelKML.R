#' @include modelApprox.R
#' @rdname interface-kml
setClass('lcModelKML', contains = 'lcApproxModel')

#. clusterTrajectories ####
#' @rdname interface-kml
#' @inheritParams clusterTrajectories
setMethod('clusterTrajectories', 'lcModelKML', function(object, at = time(object), ...) {
  if (length(at) == 0) {
    trajMat = computeKMLCenters(object)

    tsframe(
      trajMat,
      times = time(object),
      id = 'Cluster',
      time = timeVariable(object),
      response = responseVariable(object)
    )
  } else {
    callNextMethod()
  }
})


#. converged ####
#' @rdname interface-kml
setMethod('converged', 'lcModelKML', function(object) {
  TRUE
})


#' @export
#' @rdname interface-kml
logLik.lcModelKML = function(object, ...) {
  # A negated version of BIC is precomputed by kml package so let's use that
  bic = -1 * getKMLPartition(object)@criterionValues['BIC'] %>% unname()
  N = nIds(object)
  df = nClusters(object) * length(time(object)) + 1

  if (is.na(bic) && nClusters(object) == 1) {
    resp = responseVariable(object)
    # fix for kml not computing a BIC for k = 1
    alldata = cbind(model.data(object), .Mean = fitted(object)) %>%
      as.data.table()
    sigma = alldata[, sd(get(resp) - .Mean, na.rm = TRUE) * sqrt((.N - 1) / .N)]
    ll = dnorm(alldata[[resp]], alldata$.Mean, sigma, log = TRUE) %>% sum()
  }
  else {
    # recompute ll from the provided BIC
    ll = -.5 * (bic - df * log(N))
  }

  attr(ll, 'nobs') = N
  attr(ll, 'df') = df
  class(ll) = 'logLik'

  ll
}


#. postprob ####
#' @rdname interface-kml
setMethod('postprob', 'lcModelKML', function(object) {
  if (nClusters(object) == 1) {
    pp = matrix(1, nrow = nIds(object), ncol = 1)
  } else {
    ppRaw = getKMLPartition(object)@postProba
    trajClusters = kml::getClusters(object@model, nbCluster = nClusters(object), asInteger = TRUE)

    assert_that(
      length(trajClusters) == nIds(object),
      msg = sprintf('kml package returned fewer cluster assignments than expected: %d instead of %d', length(trajClusters), nIds(object))
    )

    pp = matrix(NA_real_, nrow = nIds(object), ncol = nClusters(object))
    pp[which(is.finite(trajClusters)), ] = ppRaw

    # bugfix for KML: insert missing rows for uniform postprob
    if (anyNA(pp)) {
      naMsk = is.na(pp[, 1L])
      warning(
        sprintf(
          'The kml package has outputted NA posterior probabilities for %d trajectories. Inserting uniform postprob for trajectories: \n%s',
          sum(naMsk),
          paste0('  ', ids(object)[naMsk], collapse = '\n')
        )
      )

      pp[naMsk, ] = 1 / nClusters(object)
    }
  }

  colnames(pp) = clusterNames(object)

  pp
})


#. predictPostprob
#' @rdname interface-kml
#' @inheritParams predictPostprob
setMethod('predictPostprob', 'lcModelKML', function(object, newdata, ...) {
  assert_that(
    has_name(newdata, idVariable(object)),
    has_name(newdata, timeVariable(object)),
    all(newdata[[timeVariable(object)]] %in% time(object))
  )

  valueColumn = responseVariable(object)
  centerMat = computeKMLCenters(object)
  distFun = getLcMethod(object)$distance
  times = time(object)

  if (is.null(formalArgs(distFun))) {
    affectFun = function(traj, centers) {
      kml::affectIndivC(traj = traj, clustersCenter = centers)
    }
  } else {
    affectFun = function(traj, centers) {
      kml::affectIndiv(traj = traj, clustersCenter = centers, distance = distFun)
    }
  }

  trajFun = function(trajData) {
    traj = trajData[[valueColumn]] %>% matrix(nrow=1)
    trajTimes = trajData[[timeVariable(object)]]

    affectFun(traj, centerMat[, match(trajTimes, times)]) %>%
      rep(nrow(trajData))
  }

  dtAffect = as.data.table(newdata)[, .(Cluster = trajFun(.SD)), by=c(idVariable(object))]

  pp = postprobFromAssignments(dtAffect$Cluster, nClusters(object))
  colnames(pp) = clusterNames(object)
  pp
})


getKMLPartition = function(object) {
  object@model[paste0('c', nClusters(object))][[1]]
}

computeKMLCenters = function(object) {
  trajClusters = kml::getClusters(object@model, nbCluster = nClusters(object))

  centerMat = kml::calculTrajMean(
    traj = object@model@traj,
    clust = na.omit(trajClusters),
    centerMethod = getLcMethod(object)$centerMethod
  )

  if (!is.matrix(centerMat)) {
    centerMat = matrix(centerMat, nrow = 1L)
  }
  rownames(centerMat) = clusterNames(object)

  centerMat
}
