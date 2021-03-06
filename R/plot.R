#' @include latrend.R

# plotTrajectories ####
#' @export
#' @name plotTrajectories
#' @rdname plotTrajectories
#' @aliases plotTrajectories,data.frame-method
#' @title Plot trajectories
#' @inheritParams trajectories
#' @inheritParams transformLatrendData
#' @param response Response variable `character` name or a `call`.
#' @param cluster Cluster variable name. If unspecified, trajectories are not grouped. Alternatively, cluster is a vector indicating cluster membership per id.
#' @param facet Whether to facet by cluster.
#' @examples
#' data(latrendData)
#' plotTrajectories(latrendData, response = "Y", id = "Id", time = "Time")
#'
#' plotTrajectories(latrendData, response = quote(exp(Y)), id = "Id", time = "Time")
setMethod('plotTrajectories', signature('data.frame'), function(object,
                                                                response,
                                                                time = getOption('latrend.time'),
                                                                id = getOption('latrend.id'),
                                                                cluster = NULL,
                                                                facet = TRUE,
                                                                ...) {
  if (length(cluster) > 1) {
    assert_that(length(cluster) == uniqueN(object[[id]]))
    object$Cluster = cluster[rleidv(object[[id]])]
    cluster = 'Cluster'
  }
  .plotTrajs(object, response, time, id, cluster, facet = facet)
})


.plotTrajs = function(data, response, time, id, cluster, facet = TRUE) {
  assert_that(!is.character(response) || has_name(data, response),
              has_name(data, time),
              has_name(data, id),
              is.null(cluster) || has_name(data, cluster),
              is.flag(facet))

  if (!is.null(cluster) && !facet) {
    map = aes_string(x = time, y = response, group = id, cluster = cluster, color = cluster)
  } else {
    map = aes_string(x = time, y = response, group = id, cluster = cluster)
  }

  p = ggplot(data = data, mapping = map) +
    theme(legend.position = 'top') +
    geom_line() +
    labs(title = 'Trajectories')

  if (!is.null(cluster) && facet) {
    p = p + facet_wrap(cluster)
  }
  return(p)
}

# plotClusterTrajectories ####
#' @export
#' @name plotClusterTrajectories
#' @rdname plotClusterTrajectories
#' @aliases plotClusterTrajectories,data.frame-method
#' @title Plot cluster trajectories
#' @inheritParams transformLatrendData
#' @inheritParams clusterTrajectories
#' @inheritParams predict.lcModel
#' @param object The (cluster) trajectory data.
#' @param cluster The cluster assignment column
#' @param center A function for aggregating multiple points at the same point in time
#' @param trajectories Whether to plot the original data in addition to the cluster (i.e., center) trajectories
#' @param facet Whether to facet by cluster. This is done by default when `trajectories` is enabled.
#' @param id Id column. Only needed when `trajectories = TRUE`.
#' @details Instead of passing the plotting arguments through `...`, consider modifying the ggplot2 defaults.
#' For example, changing the default line size: `update_geom_defaults("line", list(size = 1.5))`
setMethod('plotClusterTrajectories', signature('data.frame'), function(object,
    response,
    cluster = 'Cluster',
    time = getOption('latrend.time'),
    center = meanNA,
    trajectories = FALSE,
    facet = isTRUE(trajectories),
    id = getOption('latrend.id'),
    ...
  ) {
  assert_that(has_name(object, cluster),
    has_name(object, response),
    is.function(center))

  cdata = as.data.table(object) %>%
    .[, .(Value = center(get(response))), keyby=c(cluster, time)] %>%
    setnames('Value', response)

  .plotClusterTrajs(cdata,
    response = response,
    time = time,
    cluster = cluster,
    trajectories = trajectories,
    facet = facet,
    id = id,
    rawdata = object,
    ...)
})


.plotClusterTrajs = function(data, response, time, cluster = 'Cluster', trajectories = FALSE, facet = FALSE, id, rawdata = NULL, ...) {
  assert_that(
    is.data.frame(data),
    has_name(data, response),
    has_name(data, time),
    has_name(data, cluster),
    is.flag(trajectories) || is.list(trajectories),
    is.flag(facet)
  )

  if (is.factor(data[[cluster]])) {
    nClus = nlevels(data[[cluster]])
  } else {
    nClus = uniqueN(data[[cluster]])
  }

  p = ggplot(
    data = data,
    mapping = aes_string(
      x = time,
      y = response)
  )

  if (isTRUE(trajectories) || is.list(trajectories)) {
    assert_that(
      is.data.frame(rawdata),
      has_name(rawdata, id),
      has_name(rawdata, time),
      has_name(rawdata, response)
    )

    cols = c(id, time, response)
    if (facet) {
      cols = c(cols, cluster)
    }

    lineArgs = list(
      data = subset(rawdata, select = cols),
      mapping = aes_string(group = id),
      color = 'black'
    )

    if(is.list(trajectories)) {
      lineArgs = modifyList(lineArgs, trajectories)
    }

    p = p + do.call(geom_line, lineArgs)
  }

  if (facet) {
    p = p + facet_wrap(~ Cluster)
  }

  p = p + geom_line(aes_string(color = cluster), ...) +
    labs(title = 'Cluster trajectories')

  return(p)
}
