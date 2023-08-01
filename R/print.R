#' Print objects
#'
#' Print objects
#'
#' @param x an object used to select a method.
#' @param \dots further arguments passed to or from other methods.
#'
#' @name fmesher-print
NULL

#' @export
#' @param newline logical; if `TRUE` (default), end the printing with `\n`
#' @rdname fmesher-print
print.fm_segm <- function(x, newline = TRUE, ...) {
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

  cat(as.character(ret$n),
    if (isTRUE(x$is.bnd)) {
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
  if (newline) {
    cat("\n", sep = "")
  }
  return(invisible(x))
}


#' @param verbose logical
#'
#' @export
#' @rdname fmesher-print
print.fm_mesh_2d <- function(x, verbose = FALSE, ...) {
  ret <- list(verbose = verbose)
  if (verbose) {
    ret <- c(ret, list(time = x$meta$time, is.refined = x$meta$is.refined))
  }
  ret <-
    c(
      ret,
      list(
        manifold = x$manifold,
        nV = x$n,
        nT = nrow(x$graph$tv),
        xlim = if (nrow(x$loc) > 0) range(x$loc[, 1]) else NA,
        ylim = if (nrow(x$loc) > 0) range(x$loc[, 2]) else NA,
        ylim = if ((nrow(x$loc) > 0) && (ncol(x$loc) >= 3)) range(x$loc[, 3]) else NA
      )
    )
  crs <- fm_wkt(x$crs)
  ret <- c(ret, list(crs = as.character(fm_wkt(crs))))
  ret <- c(ret, list(crs_proj4 = as.character(fm_proj4string(crs))))

  if (!is.null(x$segm)) {
    ret <- c(ret, list(segm.bnd = x$segm$bnd, segm.int = x$segm$int))
  } else {
    ret <- c(ret, list(segm.bnd = fm_segm(is.bnd = TRUE), segm.int = fm_segm(is.bnd = FALSE)))
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
    print(y, ...)
    invisible(x)
  }


  cat("fm_mesh_2d object:\n", sep = "")
  cat("  CRS:\n    LegacyPROJ4:\t", ret$crs_proj4, "\n", sep = "")
  if (ret$verbose) {
    cat("    WKT:\n", ret$crs, "\n", sep = "")
  } else {
    cat("    WKT: (only shown with verbose = TRUE)", "\n", sep = "")
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
  cat("  Euler char:\t", as.character(nV - nE + nF), "\n", sep = "")

  cat("  Constraints:\t")
  print(ret$segm.bnd, newline = FALSE)
  cat(", ")
  print(ret$segm.int, newline = TRUE)
  cat("  Bounding box:\n")
  cat("    xlim: (", paste0(ret$xlim, collapse = ", "), ")", sep = "")
  cat(", ylim: (", paste0(ret$ylim, collapse = ", "), ")", sep = "")
  if (all(is.finite(ret$zlim))) {
    cat(", zlim: (", paste0(ret$zlim, collapse = ", "), ")", sep = "")
  }
  cat("\n", sep = "")
  cat("  Basis d.o.f.:\t", ret$nV, "\n", sep = "")
  invisible(x)
}



#' @param verbose logical
#'
#' @export
#' @rdname fmesher-print
print.fm_mesh_1d <- function(x, verbose = FALSE, ...) {
  cat("fm_mesh_1d object:\n", sep = "")

  cat("  Manifold:\t", x$manifold, "\n", sep = "")
  cat("  #{knots}:\t", length(x$loc), "\n", sep = "")
  cat("  Interval:\t(", paste0(x$interval, collapse = ", "), ")\n", sep = "")
  clamped <- x$free.clamped & (x$boundary == "free")
  clamped <- c("", " and clamped")[clamped + 1]
  cat("  Boundary:\t(", paste0(x$boundary, clamped, collapse = ", "), ")\n", sep = "")
  cat("  B-spline degree:\t", x$degree, "\n", sep = "")
  cat("  Basis d.o.f.:\t", x$m, "\n", sep = "")

  invisible(x)
}
