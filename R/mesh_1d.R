#' @include deprecated.R

# fm_mesh_1d ####

`match.arg.vector` <- function(arg = NULL,
                               choices,
                               length = NULL) {
  ## Like match.arg, but for a vector of options 'arg'
  if (is.null(length)) {
    length <- ifelse(is.null(arg), 1L, length(arg))
  }
  if (is.null(arg)) {
    arg <- match.arg(arg, choices)
  } else {
    for (k in seq_along(arg)) {
      arg[k] <- match.arg(arg[k], choices)
    }
  }
  if (length(arg) < length) {
    arg <- c(arg, rep(arg, length - length(arg)))
  } else if (length(arg) > length) {
    stop("Option list too long.")
  }
  return(arg)
}

#' @title Make a 1D mesh object
#' @description
#' Create a `fm_mesh_1d` object.
#'
#' @param loc B-spline knot locations.
#' @param interval Interval domain endpoints.
#' @param boundary Boundary condition specification.  Valid conditions are
#' `c('neumann', 'dirichlet', 'free', 'cyclic')`.  Two separate values can
#' be specified, one applied to each endpoint.
#' @param degree The B-spline basis degree.  Supported values are 0, 1, and 2.
#' @param free.clamped If `TRUE`, for `'free'` boundaries, clamp the
#' basis functions to the interval endpoints.
#' @param \dots Additional options, currently unused.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export
#' @family object creation and conversion
fm_mesh_1d <- function(loc,
                       interval = range(loc),
                       boundary = NULL,
                       degree = 1,
                       free.clamped = FALSE,
                       ...) {
  ## Note: do not change the order of these options without also
  ## changing 'basis.reduction' below.
  boundary.options <- c("neumann", "dirichlet", "free", "cyclic")

  boundary <- match.arg.vector(boundary, boundary.options, length = 2)
  cyclic <- !is.na(pmatch(boundary[1], "cyclic"))
  if (cyclic && is.na(pmatch(boundary[2], "cyclic"))) {
    stop("Inconsistent boundary specification 'boundary=c(",
         paste(boundary, collapse = ","), ")'.",
         sep = ""
    )
  }

  loc.orig <- loc
  if (cyclic) {
    loc <-
      (sort(unique(c(0, loc - interval[1]) %% diff(interval))) +
         interval[1])
  } else {
    loc <-
      (sort(unique(c(
        interval,
        pmax(interval[1], pmin(interval[2], loc))
      ))))
  }

  n <- length(loc)

  if (loc[1] < interval[1]) {
    stop("All 'loc' must be >= interval[1].")
  }
  if (loc[n] > interval[2]) {
    stop("All 'loc' must be <= interval[2].")
  }

  if ((degree < 0) || (degree > 2)) {
    stop(paste("'degree' must be 0, 1, or 2.  'degree=",
               degree,
               "' is not supported.",
               sep = ""
    ))
  }

  if (length(free.clamped) == 1L) {
    free.clamped <- rep(free.clamped, 2)
  }


  ## Number of basis functions
  if (degree == 0) {
    basis.reduction <- c(0, 1, 0, 1 / 2) ## neu, dir, free, cyclic
  } else if (degree == 1) {
    basis.reduction <- c(0, 1, 0, 1 / 2) ## neu, dir, free, cyclic
  } else {
    basis.reduction <- c(1, 1, 0, 1)
  }
  m <- (n + cyclic + (degree == 2) * 1
        - basis.reduction[pmatch(boundary[1], boundary.options)]
        - basis.reduction[pmatch(boundary[2], boundary.options)])
  ## if (m < 1+max(1,degree)) {
  if (m < 1L) {
    stop("Degree ", degree,
         " meshes must have at least ", 1L,
         " basis functions, not 'm=", m, "'.",
         sep = ""
    )
  }

  ## Compute representative basis midpoints.
  if ((degree == 0) || (degree == 1)) {
    mid <- loc
    if (boundary[1] == "dirichlet") {
      mid <- mid[-1]
    }
    if (boundary[2] == "dirichlet") {
      mid <- mid[-(m + 1)]
    }
  } else { ## degree==2
    if (cyclic) {
      mid <- (loc + c(loc[-1], interval[2])) / 2
    } else {
      mid <- c(loc[1], (loc[-n] + loc[-1]) / 2, loc[n])
      mid <-
        switch(boundary[1],
               neumann = mid[-1],
               dirichlet = mid[-1],
               free = mid
        )
      mid <-
        switch(boundary[2],
               neumann = mid[-(m + 1)],
               dirichlet = mid[-(m + 1)],
               free = mid
        )
    }
  }

  mesh <-
    structure(
      list(
        n = n,
        m = m,
        loc = loc,
        mid = mid,
        interval = interval,
        boundary = boundary,
        cyclic = cyclic,
        manifold = ifelse(cyclic, "S1", "R1"),
        degree = degree,
        free.clamped = free.clamped,
        idx = list(loc = NULL)
      ),
      class = c("fm_mesh_1d", "inla.mesh.1d")
    )

  if (degree < 2) {
    mesh$idx$loc <-
      fm_bary(mesh, loc.orig, method = "nearest")$index[, 1]
  } else {
    if (length(mid) >= 2) {
      mesh$idx$loc <-
        fm_bary(fm_mesh_1d(mid, degree = 0),
                loc.orig,
                method = "nearest"
        )$index[, 1]
    } else {
      mesh$idx$loc <- rep(1, length(loc.orig))
    }
  }

  return(mesh)
}



#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_mesh_1d Convert an object to `fm_mesh_1d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_mesh_1d <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_mesh_1d")
}
#' @describeIn fm_as_mesh_1d Convert each element of a list
#' @export
fm_as_mesh_1d_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "mesh_1d")
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
fm_as_mesh_1d.fm_mesh_1d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_1d", setdiff(class(x), "fm_mesh_1d"))
  x
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
#' @method fm_as_mesh_1d inla.mesh.1d
fm_as_mesh_1d.inla.mesh.1d <- function(x, ...) {
  class(x) <- c("fm_mesh_1d", class(x))
  x
}
