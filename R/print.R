#' Print objects
#'
#' Print objects
#'
#' @param x an object used to select a method.
#' @param \dots further arguments passed to or from other methods.
#' @returns The input object `x`
#' @name fmesher-print
#' @examples
#' fm_bbox(matrix(1:6, 3, 2))
#' print(fm_bbox(matrix(1:6, 3, 2)), verbose = FALSE)
#'
#' print(fmexample$mesh)
#' print(fmexample$boundary_fm)
#'
#' print(fm_mesh_1d(c(1, 2, 3, 5, 7), degree = 2))
NULL

#' @export
#' @param newline logical; if `TRUE` (default), end the printing with `\n`
#' @rdname fmesher-print
print.fm_segm <- function(x, ..., digits = NULL, verbose = TRUE, newline = TRUE) {
  my.segm <- function(x) {
    if (is.null(x)) {
      return(list(n = 0, grps = NULL))
    }
    n <- max(0, nrow(x$idx))
    if (max(0, length(unique(x$grp))) > 0) {
      grps <- unique(x$grp)
    } else {
      grps <- NULL
    }
    return(list(n = n, grps = grps))
  }

  ret <- my.segm(x)

  if (verbose) {
    cat("fm_segm object:\n", sep = "")
  }
  if (verbose && !fm_crs_is_null(fm_crs(x))) {
    crs <- fm_wkt(x$crs)
    ret <- c(ret, list(crs = as.character(fm_wkt(crs))))
    ret <- c(ret, list(crs_proj4 = as.character(fm_proj4string(crs))))

    cat("  CRS:\n    LegacyPROJ4:\t", ret$crs_proj4, "\n", sep = "")
    if (verbose) {
      if (is.na(ret$crs)) {
        cat("    WKT:\tNA\n", sep = "")
      } else {
        cat("    WKT:\n", ret$crs, "\n", sep = "")
      }
    } else {
      cat("    WKT: (only shown with verbose = TRUE)", "\n", sep = "")
    }
  }
  if (verbose) {
    cat("  ", sep = "")
  }
  cat(as.character(ret$n),
    if (all(x$is.bnd)) {
      " boundary edge"
    } else {
      " interior edge"
    },
    if (ret$n != 1) "s" else "",
    sep = ""
  )
  if (!is.null(ret$grps)) {
    n <- length(ret$grps)
    cat(" (", n, " group", if (n != 1) "s" else "", sep = "")
    if (n <= 10) {
      cat(":", ret$grps, sep = " ")
    } else {
      cat(":", ret$grps[1:10], "...", sep = " ")
    }
    cat(")")
  }
  if (verbose) {
    cat("\n  ", sep = "")
    cat("Bounding box = ", sep = "")
    print(
      fm_bbox(x),
      verbose = FALSE,
      newline = FALSE,
      digits = digits
    )
  }
  if (newline) {
    cat("\n", sep = "")
  }
  return(invisible(x))
}


#' @export
#' @param newline logical; if `TRUE` (default), end the printing with `\n`
#' @rdname fmesher-print
print.fm_segm_list <- function(x, ..., digits = NULL, verbose = FALSE, newline = TRUE) {
  if (verbose) {
    cat("list of ", length(x), " fm_segm objects:\n", sep = "")
    lapply(x, function(xx) {
      print(
        xx,
        digits = digits,
        verbose = TRUE,
        newline = TRUE
      )
    })
  } else {
    for (k in seq_along(x)) {
      print(
        x[[k]],
        digits = digits,
        verbose = FALSE,
        newline = newline
      )
      if (!newline && (k < length(x))) {
        cat(", ", sep = "")
      }
    }
  }
  return(invisible(x))
}


#' @param verbose logical
#' @param digits a positive integer indicating how many significant digits are
#' to be used for numeric and complex x. The default, NULL, uses `getOption("digits")`.
#'
#' @export
#' @rdname fmesher-print
print.fm_mesh_2d <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  if (verbose) {
    ret <- c(ret, list(time = x$meta$time, is.refined = x$meta$is.refined))
  }
  ret <-
    c(
      ret,
      list(
        manifold = x$manifold,
        nV = nrow(x$loc),
        nT = nrow(x$graph$tv)
      )
    )
  crs <- fm_wkt(x$crs)
  ret <- c(ret, list(crs = as.character(fm_wkt(crs))))
  ret <- c(ret, list(crs_proj4 = as.character(fm_proj4string(crs))))

  if (!is.null(x$segm)) {
    ret$segm <- fm_as_segm_list(list(x$segm$bnd, x$segm$int))
  } else {
    ret$segm <- fm_as_segm_list(list(
      segm.bnd = fm_segm(is.bnd = TRUE),
      segm.int = fm_segm(is.bnd = FALSE)
    ))
  }

  my.print.proc_time <- function(x, ...) {
    if (is.null(x)) {
      y <- matrix(0, 1, 5)
    } else if (!is.matrix(x)) {
      y <- matrix(x, 1, 5)
    } else {
      y <- x
    }
    for (k in seq_len(nrow(y))) {
      if (!is.na(y[k, 4L])) {
        y[k, 1L] <- y[k, 1L] + y[k, 4L]
      }
      if (!is.na(y[k, 5L])) {
        y[k, 2L] <- y[k, 2L] + y[k, 5L]
      }
    }
    y <- y[, 1L:3L, drop = FALSE]
    colnames(y) <- c(
      gettext("user"), gettext("system"),
      gettext("elapsed")
    )
    print(y, ..., digits = digits)
    invisible(x)
  }

  cat("fm_mesh_2d object:\n", sep = "")
  if (!fm_crs_is_null(x$crs)) {
    cat("  CRS:\n    LegacyPROJ4:\t", ret$crs_proj4, "\n", sep = "")
    if (ret$verbose) {
      if (is.na(ret$crs)) {
        cat("    WKT:\tNA\n", sep = "")
      } else {
        cat("    WKT:\n", ret$crs, "\n", sep = "")
      }
    } else {
      cat("    WKT: (only shown with verbose = TRUE)", "\n", sep = "")
    }
  }
  if (ret$verbose) {
    cat("  Timings:\n")
    my.print.proc_time(ret$time)
  }
  cat("  Manifold:\t", ret$manifold, "\n", sep = "")
  if (ret$verbose) {
    cat("  Refined:\t", ret$is.refined, "\n", sep = "")
  }
  nV <- ret$nV
  nE <- as.integer(sum(x$graph$vv) / 2L)
  nF <- ret$nT
  cat("  V / E / T:\t", as.character(ret$nV), " / ", sep = "")
  cat(as.character(nE), " / ", sep = "")
  cat(as.character(ret$nT), "\n", sep = "")
  cat("  Euler char.:\t", as.character(nV - nE + nF), "\n", sep = "")

  cat("  Constraints:\t")
  print(ret$segm, newline = FALSE, digits = digits, verbose = FALSE)
  cat("\n  ", sep = "")
  print(fm_bbox(x), digits = digits)
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")
  invisible(x)
}



#' @param verbose logical
#'
#' @export
#' @rdname fmesher-print
print.fm_mesh_1d <- function(x, ..., digits = NULL, verbose = FALSE) {
  cat("fm_mesh_1d object:\n", sep = "")

  cat("  Manifold:\t", x$manifold, "\n", sep = "")
  cat("  #{knots}:\t", length(x$loc), "\n", sep = "")
  cat("  Interval:\t(", paste0(format(x$interval, digits = digits),
    collapse = ", "
  ), ")\n", sep = "")
  clamped <- x$free.clamped & (x$boundary == "free")
  clamped <- c("", " and clamped")[clamped + 1]
  cat("  Boundary:\t(", paste0(x$boundary, clamped, collapse = ", "), ")\n", sep = "")
  cat("  B-spline degree:\t", x$degree, "\n", sep = "")
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")

  invisible(x)
}



#' @export
#' @rdname fmesher-print
print.fm_bbox <- function(x, ..., digits = NULL, verbose = TRUE, newline = TRUE) {
  if (verbose) {
    cat("Bounding box: ", sep = "")
  }
  if (length(x) == 0) {
    cat("NULL")
  } else {
    for (k in seq_along(x)) {
      if (k > 1) {
        cat(" x ", sep = "")
      }
      cat("(",
        paste0(format(x[[k]], digits = digits),
          collapse = ","
        ),
        ")",
        sep = ""
      )
    }
  }
  if (newline) {
    cat("\n", sep = "")
  }
  return(invisible(x))
}
