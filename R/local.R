#' @title Unit test helpers
#' @name local_testthat
#' @description Local helper functions for package unit tests
#' @param envir environment for exit handlers
#' @rdname local_testthat
#' @keywords internal
#' @returns None
#' @examples
#' outer_fun <- function() {
#'   fun <- function(envir = parent.frame()) {
#'     local_fm_testthat_assign("local_var_name", 1:4, envir = envir)
#'   }
#'   fun()
#'   local_var_name
#' }
#' exists("local_var_name")
#' outer_fun()
#' exists("local_var_name")
#'
NULL

#' @param x character; Name of variable to assign to
#' @param values the object to assign to `x`
#' @export
#' @describeIn local_testthat Assign local variable. Useful for easy cleanup
#' of global workspace with `withr::deferred_run()` when running tests
#' interactively.
local_fm_testthat_assign <- function(x, values, envir = parent.frame()) {
  exist <- exists(x, envir = envir)
  if (exist) {
    old_value <- envir[[x]]
    withr::defer(assign(x, old_value, envir = envir), envir = envir)
  } else {
    withr::defer(rm(list = x, envir = envir), envir = envir)
  }
  assign(x, values, envir = envir)
}

#' @param tolerances numeric vector of length 3; `[lowtol, midtol, hitol]`
#' @export
#' @describeIn local_testthat Assign test tolerances
#' Assign local tolerance variables. Useful for easy cleanup
#' of global workspace with `withr::deferred_run()` when running tests
#' interactively.
local_fm_testthat_tolerances <- function(tolerances = c(1e-4, 1e-2, 1e-1),
                                         envir = parent.frame()) {
  local_fm_testthat_assign("lowtol", tolerances[1], envir = envir)
  local_fm_testthat_assign("midtol", tolerances[2], envir = envir)
  local_fm_testthat_assign("hitol", tolerances[3], envir = envir)
}



#' @describeIn local_testthat Initialise environment for tests.
#' To be called either at the top of a testfile, or inside tests.
#' @export
local_fm_testthat_setup <- function(envir = parent.frame()) {
  local_fm_testthat_tolerances(envir = envir)

  sp_version <- getNamespaceVersion("sp")
  if (utils::compareVersion(sp_version, "1.6-0") >= 0) {
    if (utils::compareVersion(sp_version, "2.1-3") < 0) {
      old_sp_evolution_status <- sp::get_evolution_status()
      withr::defer(
        sp::set_evolution_status(old_sp_evolution_status),
        envir = envir
      )
    }
    fm_safe_sp(quietly = TRUE, force = TRUE)
  }

  invisible()
}









check_package_version_and_load <-
  function(pkg, minimum_version, quietly = FALSE) {
    version <- tryCatch(utils::packageVersion(pkg),
      error = function(e) NA_character_
    )
    if (is.na(version)) {
      if (!quietly) {
        message(paste0("Package '", pkg, "' is not installed."))
      }
      return(NA_character_)
    }
    if (version < minimum_version) {
      if (!quietly) {
        message(paste0(
          "Installed '", pkg, "' version is ", version, " but ",
          "version >= ", minimum_version, " is required."
        ))
      }
      return(NA_character_)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!quietly) {
        message("Package '", pkg, "' not loaded safely.")
      }
      return(NA_character_)
    }
    return(version)
  }


#' Check for potential `sp` version compatibility issues
#'
#' Loads the sp package with `requireNamespace("sp", quietly = TRUE)`, and
#' checks and optionally sets the `sp` evolution status flag if `rgdal` is unavailable.
#' This function is only needed for backwards compatibility with `sp` versions
#' before `2.0-0`.
#'
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default `FALSE`
#' @param force logical; If `rgdal` is unavailable
#' and evolution status is less that `2L`, return `FALSE` if `force` is `FALSE`.
#' If `force` is `TRUE`, return `TRUE` if the package configuration is safe,
#' potentially after forcing the evolution status to `2L`.
#' Default `FALSE`
#' @param minimum_version character; the minimum required sp version.
#' Default 1.4-5 (should always match the requirement in the package
#' DESCRIPTION)
#' @return Returns (invisibly) `FALSE` if a potential issue is detected, and give a
#' message if `quietly` is `FALSE`. Otherwise returns `TRUE`
#' @export
#' @examples
#' if (fm_safe_sp()) {
#'   # Run sp dependent calculations
#' }
#' @keywords internal
fm_safe_sp <- function(quietly = FALSE,
                       force = FALSE,
                       minimum_version = "1.4-5") {
  sp_version <-
    check_package_version_and_load(
      pkg = "sp",
      minimum_version = minimum_version,
      quietly = quietly
    )
  if (is.na(sp_version)) {
    return(invisible(FALSE))
  }

  if (sp_version >= "2.1.4") {
    return(invisible(TRUE))
  }

  if (sp_version >= "1.6-0") {
    # Default to 2L to allow future sp to stop supporting
    # get_evolution_status; assume everything is fine if it fails.
    if (sp_version < "2.1-3") {
      # From at least version 2.1-3, sp gives a warning when using
      # get_evolution_status, and the status is fixed to 2L, so no
      # need to set it.
      evolution_status <- tryCatch(
        sp::get_evolution_status(),
        error = function(e) 2L
      )
    } else {
      evolution_status <- 2L
    }
    rgdal_version <- tryCatch(utils::packageVersion("rgdal"),
      error = function(e) NA_character_
    )
    if ((evolution_status < 2L) && is.na(rgdal_version)) {
      if (!quietly) {
        message("'sp' version >= 1.6-0 detected, rgdal isn't installed, and evolution status is < 2L.")
      }
      if (!force) {
        if (!quietly) {
          message(
            "This may cause issues with some CRS handling code. To avoid this, use 'sp::set_evolution_status(2L)'"
          )
        }
        return(invisible(FALSE))
      }

      sp::set_evolution_status(2L)
      if (!quietly) {
        message(
          "Ran 'sp::set_evolution_status(2L)' to avoid issues with some CRS handling code."
        )
      }
    }
  }
  return(invisible(TRUE))
}
