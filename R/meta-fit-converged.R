#' @include meta-fit.R

#' @export
#' @rdname lcFitMethods
#' @examples
#'
#' data(latrendData)
#' method <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time", nClusters = 2)
#' metaMethod <- lcFitConverged(method, maxRep = 10)
#' metaMethod
#' model <- latrend(metaMethod, latrendData)
setClass('lcFitConverged', contains = 'lcMetaMethod')

#' @export
#' @rdname lcFitMethods
#' @param method The `lcMethod` to use for fitting.
#' @param maxRep The maximum number of fit attempts
lcFitConverged = function(method, maxRep = Inf) {
  mc = match.call.all()
  mc$method = NULL
  mc$Class = 'lcFitConverged'

  object = do.call(new, as.list(mc))
  object@method = method
  object
}


#' @rdname interface-metaMethods
setMethod('fit', 'lcFitConverged', function(method, data, envir, verbose) {
  attempt = 1L

  repeat {
    enter(verbose, level = verboseLevels$fine, suffix = '')
    model = fit(getLcMethod(method), data = data, envir = envir, verbose = verbose)
    exit(verbose, level = verboseLevels$fine, suffix = '')

    if (converged(model)) {
      return (model)
    } else if (attempt >= method$maxRep) {
      warning(
        sprintf(
          'Failed to obtain converged result (got %s) for %s within %d attempts.\n\tReturning last model.',
          converged(model),
          class(getLcMethod(method))[1],
          method$maxRep
        ),
        immediate. = TRUE
      )
      return (model)
    } else {
      attempt = attempt + 1L

      if (has_lcMethod_args(getLcMethod(method), 'seed')) {
        seed = sample.int(.Machine$integer.max, 1L)
        set.seed(seed)
        # update fit method with new seed
        method@method = update(getLcMethod(method), seed = seed, .eval = TRUE)
      }

      if (is.infinite(method$maxRep)) {
        cat(
          verbose,
          sprintf(
            'Method failed to converge (got %gs). Retrying... attempt %d',
            converged(model),
            attempt
          )
        )
      } else {
        cat(
          verbose,
          sprintf(
            'Method failed to converge (got %s). Retrying... attempt %d / %d',
            converged(model),
            attempt,
            method$maxRep
          )
        )
      }
    }
  }
})

#' @rdname interface-metaMethods
setMethod('validate', 'lcFitConverged', function(method, data, envir = NULL, ...) {
  callNextMethod()

  validate_that(
    has_lcMethod_args(method, 'maxRep'),
    is.count(method$maxRep)
  )
})
