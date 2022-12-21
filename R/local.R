#' @title Unit test helpers
#' @name local_testthat
#' @description Local helper functions for package unit tests
#' @param envir environment for exit handlers
#' @rdname local_testthat
#' @keywords internal
NULL

#' @param x character; Name of variable to assign to
#' @param values the object to assign to `x`
#' @export
#' @describeIn local_testthat Assign local variable. Useful for easy cleanup
#' of global workspace with `withr::deferred_run()` when running tests
#' interactively.
local_testthat_assign <- function(x, values, envir = parent.frame()) {
  exist <- exists(x, envir = envir)
  if (exist) {
    old_value <- envir[[x]]
    withr::defer(assign(x, old_value, envir = envir), envir = envir)
  } else {
    withr::defer(rm(list = x, envir = envir), envir = envir)
  }
  assign(x, values, envir = envir)
}






#' @describeIn local_testthat Tests should set num.threads = "1:1" to ensure
#' within-system repeatability by calling `local_fm_safe_inla()`;
#' see also [fm_safe_inla()]
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered. Default: `FALSE`, multicore
#' not allowed (used for examples and unit tests).
#' @param quietly logical; if `TRUE`, prints diagnostic messages. A message is
#' always printed if the INLA `num.threads` option is altered, regardless of the
#' `quietly` argument. Default: TRUE.
#' @export
local_fm_safe_inla <- function(multicore = FALSE,
                               quietly = TRUE,
                               envir = parent.frame()) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    # Save the num.threads option so it can be restored
    old_threads <- tryCatch(
      INLA::inla.getOption("num.threads"),
      error = function(e) {
        e
      }
    )
    if (inherits(old_threads, "simpleError")) {
      return(testthat::skip("inla.getOption() failed, skip INLA tests."))
    }
    withr::defer(
      INLA::inla.setOption(num.threads = old_threads),
      envir
    )
    # Save the fmesher.timeout option so it can be restored
    old_fmesher_timeout <- INLA::inla.getOption("fmesher.timeout")
    withr::defer(
      INLA::inla.setOption(fmesher.timeout = old_fmesher_timeout),
      envir
    )
    INLA::inla.setOption(fmesher.timeout = 30)
  }
  testthat::skip_if_not(fm_safe_inla(multicore = multicore, quietly = quietly))
}


#' @describeIn local_testthat Initialise environment for tests.
#' To be called either at the top of a testfile, or inside tests.
#' Does *not* call [local_fm_safe_inla()], since that may invoke a skip and
#' should be called inside each test that relies on INLA.
#' @export
local_fm_testthat_setup <- function(envir = parent.frame()) {
  # Nothing to do yet
}







#' Load INLA safely for examples and tests
#'
#' Loads the INLA package with `requireNamespace("INLA", quietly = TRUE)`, and
#' optionally checks and sets the multicore `num.threads` INLA option.
#'
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered.
#' If `FALSE`, forces `num.threads="1:1"`. Default: NULL, checks
#' if running in testthat or non-interactively, in which case sets
#' `multicore=FALSE`, otherwise `TRUE`.
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default: FALSE.
#' @export
#' @return logical; `TRUE` if INLA was loaded safely, otherwise FALSE
#'
#' @examples
#' \dontrun{
#' if (fm_safe_inla()) {
#'   # Run inla dependent calculations
#' }
#' }
#'
fm_safe_inla <- function(multicore = NULL,
                         quietly = FALSE) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    if (is.null(multicore)) {
      multicore <-
        !identical(Sys.getenv("TESTTHAT"), "true") ||
          interactive()
    }
    if (!multicore) {
      n.t <- INLA::inla.getOption("num.threads")
      if (!quietly) {
        message(paste0("Current num.threads is '", n.t, "'."))
      }
      if (!identical(n.t, "1:1")) {
        if (!quietly) {
          message(paste0(
            "Setting INLA option num.threads to '1:1'.",
            " Previous value '", n.t, "'."
          ))
        }
        INLA::inla.setOption(num.threads = "1:1")
      } else {
        if (!quietly) {
          message("No num.threads change needed.")
        }
      }
    }
    TRUE
  } else {
    if (!quietly) {
      message("INLA not loaded safely.")
    }
    FALSE
  }
}
