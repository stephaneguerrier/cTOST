#' Restore the caller's random number stream on exit
#'
#' Several functions fix the seed internally so that their Monte Carlo results
#' are reproducible. Calling this helper at the top of such a function records
#' the current state of R's random number generator and restores it when the
#' function exits, so that the caller's own stream is not altered.
#'
#' @param envir The frame whose exit triggers the restore (default: the caller).
#' @return Invisibly, \code{NULL}.
#' @keywords internal
#' @noRd
.restore_rng_on_exit = function(envir = parent.frame()) {
  has_seed = exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed = if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
  expr = if (has_seed) {
    bquote(assign(".Random.seed", .(old_seed), envir = globalenv()))
  } else {
    quote(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    })
  }
  do.call(on.exit, list(expr, add = TRUE), envir = envir)
  invisible(NULL)
}

#' Check that an argument is a single finite number within bounds
#'
#' @param x The value to check.
#' @param name The argument name used in the error message.
#' @param min,max Allowed range.
#' @param strict_min,strict_max Whether the bounds are excluded.
#' @param integer Whether the value must be integer-valued.
#' @return Invisibly, \code{x}.
#' @keywords internal
#' @noRd
.check_scalar = function(x, name, min = -Inf, max = Inf, strict_min = FALSE,
                         strict_max = FALSE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite number.", name), call. = FALSE)
  }
  if (x < min || (strict_min && x == min)) {
    stop(sprintf("'%s' must be %s %s.", name, if (strict_min) ">" else ">=", format(min)), call. = FALSE)
  }
  if (x > max || (strict_max && x == max)) {
    stop(sprintf("'%s' must be %s %s.", name, if (strict_max) "<" else "<=", format(max)), call. = FALSE)
  }
  if (integer && x != round(x)) {
    stop(sprintf("'%s' must be an integer.", name), call. = FALSE)
  }
  invisible(x)
}

#' Resolve a seed argument for the differential-privacy tests
#'
#' @param seed \code{NULL} (draw a seed from the caller's stream) or a single
#'   integer-valued number.
#' @param margin Reserve so that \code{seed + margin} is still an integer.
#' @return An integer seed.
#' @keywords internal
#' @noRd
.resolve_seed = function(seed, margin = 1000L) {
  if (is.null(seed)) {
    return(sample.int(.Machine$integer.max - margin, 1))
  }
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed != round(seed) ||
      abs(seed) > .Machine$integer.max - margin) {
    stop(sprintf("'seed' must be NULL or a single integer with absolute value at most %d.",
                 .Machine$integer.max - margin), call. = FALSE)
  }
  as.integer(seed)
}
