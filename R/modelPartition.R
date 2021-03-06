#' @include modelApprox.R
setClass(
  'lcModelPartition',
  representation(
    center = 'function',
    clusterTrajectories = 'data.frame',
    postprob = 'matrix',
    name = 'character'
  ),
  contains = 'lcApproxModel'
)

#' @export
#' @title Create a lcModel with pre-defined partitioning
#' @description Represents an arbitrary partitioning of a set of trajectories.
#' As such, this model has no predictive capabilities. The cluster trajectories are represented by the specified center function (mean by default).
#' @inheritParams lcMethodStratify
#' @param data A `data.frame` representing the trajectory data.
#' @param trajectoryAssignments A `vector` of cluster membership per trajectory, either `factor`, or `integer` (`1` to `nClusters`).
#' @param nClusters The number of clusters. Optional for `factor` assignments.
#' @param clusterNames The names of the clusters, or a function with input `n` outputting a `character vector` of names.
#' @param envir The `environment` associated with the model. Used for evaluating the assigned `data` object by [model.data.lcModel].
lcModelPartition = function(data,
                            response,
                            trajectoryAssignments,
                            nClusters = NA,
                            center = meanNA,
                            clusterNames = NULL,
                            time = getOption('latrend.time'),
                            id = getOption('latrend.id'),
                            name = 'part',
                            envir = parent.frame()) {
  assert_that(
    is.data.frame(data),
    has_name(data, response),
    has_name(data, time),
    has_name(data, id),
    is.character(clusterNames) || is.null(clusterNames),
    length(clusterNames) %in% c(0, nClusters),
    is.function(center)
  )

  assert_that(
    all(vapply(
      trajectoryAssignments, is.count, FUN.VALUE = TRUE
    )) || is.factor(trajectoryAssignments),
    length(trajectoryAssignments) == uniqueN(data[[id]])
  )

  if (is.factor(trajectoryAssignments)) {
    assert_that(is.na(nClusters) || nlevels(trajectoryAssignments) == nClusters)
  }

  intAssignments = as.integer(trajectoryAssignments)
  assert_that(is.na(nClusters) || max(intAssignments) <= nClusters)

  # Determine number of clusters
  if (is.na(nClusters)) {
    if (is.factor(trajectoryAssignments)) {
      numClus = nlevels(trajectoryAssignments)
    } else {
      numClus = max(intAssignments)
    }
  } else {
    numClus = nClusters
  }
  assert_that(min(intAssignments) >= 1, max(intAssignments) <= numClus)

  pp = postprobFromAssignments(intAssignments, k = numClus)

  if (is.null(clusterNames)) {
    clusterNames = make.clusterNames(numClus)
  }

  clusTrajs = computeCenterClusterTrajectories(
    data,
    assignments = intAssignments,
    nClusters = numClus,
    fun = center,
    id = id,
    time = time,
    response = response
  )

  mc = match.call()
  model = new(
    'lcModelPartition',
    call = mc,
    data = data,
    center = center,
    clusterTrajectories = clusTrajs,
    postprob = pp,
    name = name,
    clusterNames = clusterNames,
    id = id,
    time = time,
    response = response
  )
  environment(model) = envir
  return(model)
}

#' @rdname interface-custom
setMethod('clusterTrajectories', signature('lcModelPartition'), function(object, at = time(object), ...) {
  if (is.null(at)) {
    clusTrajs = as.data.table(object@clusterTrajectories)
    clusTrajs[, Cluster := factor(Cluster,
                                  levels = seq_len(nClusters(object)),
                                  labels = clusterNames(object))]
    return(clusTrajs[])
  } else {
    callNextMethod()
  }
})


#. converged ####
#' @rdname interface-custom
setMethod('converged', signature('lcModelPartition'), function(object, ...) {
  TRUE
})


# . getName ####
#' @rdname interface-custom
setMethod('getName', signature('lcModelPartition'), function(object, ...) object@name)

# . getShortName ####
#' @rdname interface-custom
setMethod('getShortName', signature('lcModelPartition'), function(object, ...)
  object@name)


#. postprob ####
#' @rdname interface-custom
setMethod('postprob', signature('lcModelPartition'), function(object, ...) {
  pp = object@postprob
  colnames(pp) = clusterNames(object)
  return(pp)
})



computeCenterClusterTrajectories = function(data,
                                            assignments,
                                            nClusters,
                                            fun = mean,
                                            id,
                                            time,
                                            response) {
  assert_that(
    is.data.frame(data),
    has_name(data, response),
    has_name(data, time),
    has_name(data, id)
  )
  assert_that(nClusters >= 1)
  assert_that(
    is.integer(assignments),
    all(is.finite(assignments)),
    all(vapply(assignments, is.count, FUN.VALUE = TRUE)),
    length(assignments) == uniqueN(data[[id]]),
    min(assignments) >= 1,
    max(assignments) <= nClusters
  )
  assert_that(is.function(fun))

  rowClusters = assignments[rleidv(data[[id]])]
  clusTrajs = as.data.table(data) %>%
    .[, .(Value = fun(get(response))), by = .(Cluster = rowClusters, Time = get(time))]

  if (uniqueN(assignments) < nClusters) {
    warning(
      'empty clusters present. cluster trajectory for empty clusters will be set constant at 0'
    )
    # add missing clusters
    emptyClusTraj = clusTrajs[, .(Time = unique(Time), Value = 0)]
    clusTrajs = rbind(clusTrajs,
                      data.table(Cluster = rep(
                        setdiff(seq_len(nClusters), unique(rowClusters)), each = nrow(emptyClusTraj)
                      ),
                      emptyClusTraj))
  }

  setnames(clusTrajs, 'Value', response)
  setnames(clusTrajs, 'Time', time)
  return(clusTrajs[])
}
