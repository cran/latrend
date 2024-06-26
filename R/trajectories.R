#' @include generics.R latrend.R


# trajectories ####
#' @rdname trajectories
#' @param cluster Experimental feature for data.frame input: a vector of cluster membership per id
#' @aliases trajectories,data.frame-method
setMethod('trajectories', 'data.frame', function(object, id, time, response, cluster, ...) {
  if (length(cluster) > 1 && length(cluster) == uniqueN(object[[id]])) {
    object$Cluster = cluster[rleidv(object[[id]])]
  }

  object
})


#' @rdname trajectories
#' @aliases trajectories,matrix-method
setMethod('trajectories', 'matrix', function(object, id, time, response, ...) {
  data = tsframe(object, id = id, time = time, response = response)
  trajectories(data, id = id, time = time, response = response, ...)
})


#' @rdname trajectories
#' @aliases trajectories,call-method
#' @param envir The `environment` used to evaluate the data object in (e.g., in case `object` is of type `call`).
setMethod('trajectories', 'call', function(object, ..., envir) {
  data = eval(object, envir = envir)
  trajectories(data, ..., envir = envir)
})



# plotCusterTrajectories ####
#' @export
#' @name plotClusterTrajectories
#' @rdname plotClusterTrajectories
#' @aliases plotClusterTrajectories,data.frame-method
#' @title Plot cluster trajectories
#' @inheritParams trajectories
#' @inheritParams clusterTrajectories
#' @inheritParams predict.lcModel
#' @param object The (cluster) trajectory data.
#' @param cluster The cluster assignment column
#' @param center A function for aggregating multiple points at the same point in time
#' @param trajectories Whether to additionally plot the original trajectories (`TRUE`),
#' or to show the expected interval (standard deviation, standard error, range, or percentile range)
#' of the observations at the respective moment in time.
#'
#' Note that visualizing the expected intervals is currently only supported for time-aligned trajectories,
#' as the interval is computed at each unique moment in time.
#' By default (`FALSE`), no information on the underlying trajectories is shown.
#' @param facet Whether to facet by cluster. This is done by default when `trajectories` is enabled.
#' @param id Id column. Only needed when `trajectories = TRUE`.
setMethod('plotClusterTrajectories', 'data.frame', function(
  object,
  response,
  cluster = 'Cluster',
  clusterOrder = character(),
  clusterLabeler = make.clusterPropLabels,
  time = getOption('latrend.time'),
  center = meanNA,
  trajectories = c(FALSE, 'sd', 'se', '80pct', '90pct', '95pct', 'range'),
  facet = !isFALSE(as.logical(trajectories[1])),
  id = getOption('latrend.id'),
  ...
) {
  assert_that(
    has_name(object, cluster),
    has_name(object, response),
    is.function(center)
  )

  clusTrajData = as.data.table(object) %>%
    .[, .(Value = center(get(..response))), keyby = c(cluster, time)] %>%
    setnames('Value', response)

  clusterNames = as.character(unique(clusTrajData[[cluster]]))
  clusterSizes = clusTrajData[, uniqueN(id), by = c(cluster)]$V1
  clusterLabels = clusterLabeler(clusterNames, clusterSizes)


  .plotClusterTrajs(
    clusTrajData,
    clusterOrder = clusterOrder,
    clusterLabels = clusterLabels,
    response = response,
    time = time,
    cluster = cluster,
    trajectories = trajectories[1],
    facet = facet,
    id = id,
    trajData = object,
    ...
  )
})


.plotClusterTrajs = function(
  data,
  clusterOrder,
  clusterLabels,
  response,
  time,
  cluster = 'Cluster',
  trajectories = FALSE,
  facet = FALSE,
  id,
  trajData = NULL,
  ...
) {
  .loadOptionalPackage('ggplot2')

  assert_that(
    is.data.frame(data),
    has_name(data, response),
    has_name(data, time),
    has_name(data, cluster),
    is.scalar(trajectories),
    is.flag(trajectories) || is.character(trajectories),
    is.flag(facet)
  )

  clusterNames = as.character(unique(data[[cluster]]))
  clusterOrderNames = make.orderedClusterNames(clusterNames, clusterOrder, subset = TRUE)
  clusterIndices = match(clusterOrderNames, clusterNames)

  data = data.table(data)[get(cluster) %in% clusterOrderNames] %>%
    .[, c(cluster) := factor(
      get(cluster),
      levels = clusterNames[clusterIndices],
      labels = clusterLabels[clusterIndices]
    )]

  if (!is.null(trajData)) {
    trajData = data.table(trajData)[get(cluster) %in% clusterOrderNames] %>%
      .[, c(cluster) := factor(
        get(cluster),
        levels = clusterNames[clusterIndices],
        labels = clusterLabels[clusterIndices]
      )]
  }


  if (isTRUE(as.logical(trajectories))) {
    # show trajectories
    assert_that(
      is.data.frame(trajData),
      nrow(trajData) > 0,
      !missing(id)
    )

    p = plotTrajectories(
      trajData,
      response = response,
      time = time,
      id = id,
      facet = FALSE,
      cluster = cluster,
      ...
    )
  } else if (isFALSE(as.logical(trajectories))) {
    # don't show trajectories
    p = ggplot2::ggplot()
  } else {
    # ribbon plot
    p = ggplot2::ggplot()

    if (grepl('^[0-9]{1,2}pct$', trajectories)) {
      # percentile range
      pct = as.integer(substr(trajectories, 0, nchar(trajectories) - 3))
      assert_that(
        pct > 0,
        msg = 'invalid plotClusterTrajectories() argument trajectories = "##pct": ## must be > 0'
      )

      ribbonFun = function(x) quantile(x, probs = .5 + c(-1, 1) * pct / 200, na.rm = TRUE)
      ribbonTitle = sprintf('%dth percentile interval', pct)
    } else {
      trajectories = match.arg(
        trajectories,
        choices = c('sd', 'se', '80pct', '90pct', '95pct', 'range')
      )

      se = function(x) {
        sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
      }

      ribbonFun = switch(trajectories,
        sd = function(x) c(-1, 1) * sd(x, na.rm = TRUE) + mean(x, na.rm = TRUE),
        se = function(x) c(-1, 1) * se(x) + mean(x, na.rm = TRUE),
        range = function(x) range(x, na.rm = TRUE)
      )

      ribbonTitle = switch(trajectories,
        sd = 'SD',
        se = 'standard error',
        range = 'range'
      )
    }

    # compute ribbon
    ribbonData = as.data.table(trajData)[, as.list(ribbonFun(get(response))), keyby = c(cluster, time)]
    assert_that(
      ncol(ribbonData) == 4,
      msg = 'ribbon stat implementation error\nplease report.'
    )
    setnames(ribbonData, c(cluster, time, 'ymin', 'ymax'))

    p = p + ggplot2::geom_ribbon(
      mapping = ggplot2::aes(
        x = !!as.name(time),
        ymin = ymin,
        ymax = ymax,
        fill = !!.as_lang(cluster)
      ),
      data = ribbonData,
      alpha = .5
    ) + ggplot2::labs(
      subtitle = sprintf(
        'Range represents %s of local trajectory observations',
        ribbonTitle
      )
    )
  }

  if (facet) {
    # black cluster trajectories
    p = p +
      ggplot2::guides(color = 'none') +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = !!.as_lang(time),
          y = !!.as_lang(response),
          group = !!.as_lang(cluster)
        ),
        data = data,
        ...
      )
  } else {
    # non-faceted plot
    if (trajectories) {
      # color used for trajectories, use shapes to mark cluster trajectories
      p = p + ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = !!.as_lang(time),
          y = !!.as_lang(response),
          group = !!.as_lang(cluster)
        ),
        data = data,
        ...
      )
    } else {
      # colored cluster trajectories
      p = p + ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = !!.as_lang(time),
          y = !!.as_lang(response),
          color = !!.as_lang(cluster)
        ),
        data = data,
        ...
      )
    }
  }

  if (facet) {
    p = p + ggplot2::facet_wrap(cluster)
  }

  p + ggplot2::labs(title = 'Cluster trajectories')
}


# plotTrajectories ####
#' @export
#' @name plotTrajectories
#' @rdname plotTrajectories
#' @aliases plotTrajectories,data.frame-method
#' @inheritParams trajectories
#' @param response Response variable `character` name or a `call`.
#' @param cluster Whether to plot trajectories grouped by cluster (determined by the "Cluster" column).
#' Alternatively, the name of the cluster column indicating trajectory cluster membership.
#' If unspecified, trajectories are grouped if the object contains a "Cluster" column.
#' @param facet Whether to facet by cluster.
#' @seealso [trajectories] [plotFittedTrajectories] [plotClusterTrajectories]
#' @examples
#' data(latrendData)
#'
#' if (require("ggplot2")) {
#'   plotTrajectories(latrendData, response = "Y", id = "Id", time = "Time")
#'
#'   plotTrajectories(
#'     latrendData,
#'     response = quote(exp(Y)),
#'     id = "Id",
#'     time = "Time"
#'   )
#'
#'   plotTrajectories(
#'     latrendData,
#'     response = "Y",
#'     id = "Id",
#'     time = "Time",
#'     cluster = "Class"
#'   )
#' }
setMethod('plotTrajectories', 'data.frame',
  function(
    object,
    response,
    cluster,
    time = getOption('latrend.time'),
    id = getOption('latrend.id'),
    facet = TRUE,
    ...
  ) {
  .loadOptionalPackage('ggplot2')

  assert_that(
    has_name(object, time),
    has_name(object, id)
  )

  if (missing(cluster)) {
    cluster = has_name(object, 'Cluster')
  }

  if (isTRUE(cluster)) {
    assert_that(
      has_name(object, 'Cluster'),
      msg = 'cluster = TRUE but object has no "Cluster" column'
    )
    cluster = 'Cluster'
  }

  if (missing(response)) {
    response = .guessResponseVariable(object, id = id, time = time, cluster = cluster)
  }

  assert_that(
    !is.character(response) || has_name(object, response),
    is.scalar(cluster),
    isFALSE(cluster) || is.character(cluster),
    is.flag(facet)
  )

  if (is.character(cluster)) {
    assert_that(has_name(object, cluster))
  }

  if (!isFALSE(cluster) && !isTRUE(facet)) {
    map = ggplot2::aes(
      x = !!.as_lang(time),
      y = !!.as_lang(response),
      group = !!.as_lang(id),
      color = !!.as_lang(cluster)
    )
  } else {
    map = ggplot2::aes(
      x = !!.as_lang(time),
      y = !!.as_lang(response),
      group = !!.as_lang(id)
    )
  }

  p = ggplot2::ggplot() +
    ggplot2::theme(legend.position = 'top') +
    ggplot2::geom_line(mapping = map, data = object) +
    ggplot2::labs(title = 'Trajectories')

  if (!isFALSE(cluster) && isTRUE(facet)) {
    p = p + ggplot2::facet_wrap(cluster)
  }

  p
})


#' @export
#' @rdname plotTrajectories
#' @aliases plotTrajectories,ANY-method
setMethod('plotTrajectories', 'ANY', function(object, ...) {
  data = trajectories(object, ...)
  plotTrajectories(data, ...)
})
