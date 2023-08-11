#' @title Draw `fm_segm` objects.
#'
#' @description
#' Draws a [fm_segm()] object with generic or `rgl`
#' graphics.
#'
#' @importFrom grDevices cm.colors rgb
#' @importFrom graphics lines plot.window points
#'
#' @param x An [fm_segm()] object.
#' @param loc Point locations to be used if `x$loc` is `NULL`.
#' @param col Segment color specification.
#' @param colors Colors to cycle through if `col` is `NULL`.
#' @param add If `TRUE`, add to the current plot, otherwise start a new
#' plot.
#' @param xlim,ylim X and Y axis limits for a new plot.
#' @param rgl If `TRUE`, use `rgl` for plotting.
#' @param asp Aspect ratio for new plots. Default 1.
#' @param axes logical; whether axes should be drawn on the plot.
#' Default FALSE.
#' @param xlab,ylab character; labels for the axes.
#' @param \dots Additional parameters, passed on to graphics methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_segm()]
#' @export
plot.fm_segm <- function(x, ..., add = FALSE) {
  lines(x, add = add, rgl = FALSE, ...)
}

#' @rdname plot.fm_segm
#' @export
lines.fm_segm <- function(x, loc = NULL, col = NULL,
                          colors = c("black", "blue", "red", "green"),
                          add = TRUE, xlim = NULL, ylim = NULL,
                          rgl = FALSE, asp = 1,
                          axes = FALSE,
                          xlab = "",
                          ylab = "",
                          ...) {
  segm <- x
  if (!is.null(segm$loc)) {
    loc <- segm$loc
  }
  stopifnot(!is.null(loc), ncol(loc) >= 2)
  if (ncol(loc) < 3) {
    loc <- cbind(loc, 0.0)
  }
  color <- col
  dev <- NULL
  if (rgl) {
    requireNamespace("rgl", quietly = FALSE)
    if (!add) {
      dev <- rgl::open3d()
      rgl::view3d(0, 0, fov = 0)
    } else {
      dev <- rgl::cur3d()
    }
  } else if (!add) {
    idx <- unique(as.vector(segm$idx))
    if (is.null(xlim)) {
      xlim <- range(loc[idx, 1])
    }
    if (is.null(ylim)) {
      ylim <- range(loc[idx, 2])
    }
    plot(NA, type = "n",
         xlim = xlim, ylim = ylim, asp = asp,
         axes = axes,
         xlab = xlab, ylab = ylab,
         ...)
  }

  grps <- if (is.null(segm$grp)) rep(0L, nrow(segm$idx)) else segm$grp
  for (grp in unique(grps)) {
    idx <- which(grps == grp)
    if (is.null(col)) {
      color <- colors[1 + (grp %% length(colors))]
    }
    if (rgl) {
      rgl::segments3d(loc[as.vector(t(segm$idx[idx, , drop = FALSE])), , drop = FALSE],
        color = color,
        ...
      )
    } else {
      lines(loc[t(cbind(segm$idx[idx, , drop = FALSE], NA)), 1],
        loc[t(cbind(segm$idx[idx, , drop = FALSE], NA)), 2],
        col = color,
        ...
      )
    }
  }
  return(invisible(dev))
}

#' Generate text RGB color specifications.
#'
#' Generates a tex RGB color specification matrix based on a color palette.
#'
#' @keywords internal
#' @param color `character`, `matrix` or `vector`
#' @param color.axis The min/max limit values for the color mapping.
#' @param color.n The number of colors to use in the color palette.
#' @param color.palette A color palette function.
#' @param color.truncate If `TRUE`, truncate the colors at the color axis
#' limits.
#' @param alpha Transparency/opaqueness values.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
generate_colors <- function(color,
                            color.axis = NULL,
                            color.n = 512,
                            color.palette = cm.colors,
                            color.truncate = FALSE,
                            alpha = NULL) {
  if (is.character(color)) {
    colors <- color
  } else if (is.vector(color) || (is.matrix(color) && (ncol(color) == 1))) {
    if (is.null(color.axis)) {
      color.axis <- c(min(color, na.rm = TRUE), max(color, na.rm = TRUE))
    }
    if (color.truncate) {
      not.ok <- ((color < color.axis[1]) |
        (color > color.axis[2]))
    } else {
      not.ok <- rep(FALSE, length(color))
    }
    cs <- (pmax(color.axis[1],
      pmin(color.axis[2], color, na.rm = TRUE),
      na.rm = TRUE
    ))
    cs <- (cs - color.axis[1]) / (color.axis[2] - color.axis[1])
    not.ok <- not.ok | is.na(cs)
    cs[not.ok] <- 0.5
    if (is.null(alpha)) {
      alpha <- as.numeric(!not.ok)
    } else {
      alpha[not.ok] <- 0
    }

    ics <- (as.numeric(cut(cs, seq(0, 1, length.out = color.n + 1),
      include.lowest = TRUE
    )))
    colors <- color.palette(color.n)[ics]

    ## Todo: handle alpha, combining "input alpha" with "not.ok-alpha"
  } else if (is.matrix(color) && (ncol(color) == 3)) {
    if (is.null(color.axis)) {
      color.axis <- c(min(color, na.rm = TRUE), max(color, na.rm = TRUE))
    }
    if (color.truncate) {
      not.ok <- ((color[, 1] < color.axis[1]) |
        (color[, 2] < color.axis[1]) |
        (color[, 3] < color.axis[1]) |
        (color[, 1] > color.axis[2]) |
        (color[, 2] > color.axis[2]) |
        (color[, 3] > color.axis[2]))
    } else {
      not.ok <- rep(FALSE, nrow(color))
    }
    cs <- matrix(
      pmax(color.axis[1],
        pmin(color.axis[2], color, na.rm = TRUE),
        na.rm = TRUE
      ), dim(color)
    )
    cs <- (cs - color.axis[1]) / (color.axis[2] - color.axis[1])
    not.ok <- not.ok | is.na(cs[, 1]) | is.na(cs[, 2]) | is.na(cs[, 3])
    cs[not.ok, ] <- c(0.5, 0.5, 0.5)
    if (is.null(alpha)) {
      alpha <- as.numeric(!not.ok)
    } else {
      alpha[not.ok] <- 0
    }
    colors <- rgb(cs[, 1], cs[, 2], cs[, 3])
  } else {
    stop("color specification must be character, matrix, or vector.")
  }

  return(list(colors = colors, alpha = alpha))
}



#' Low level triangulation mesh plotting
#'
#' Plots a triangulation mesh using `rgl`.
#'
#'
#' @param x A `fm_mesh_2d()` object
#' @param col Color specification.  A single named color, a vector of scalar
#' values, or a matrix of RGB values.
#' @param color.axis The min/max limit values for the color mapping.
#' @param color.n The number of colors to use in the color palette.
#' @param color.palette A color palette function.
#' @param color.truncate If `TRUE`, truncate the colors at the color axis
#' limits.
#' @param alpha Transparency/opaqueness values. See `rgl.material`.
#' @param lwd Line width for edges. See `rgl.material`.
#' @param specular Specular color. See `rgl.material`.
#' @param draw.vertices If `TRUE`, draw triangle vertices.
#' @param size Size for vertex points.
#' @param draw.edges If `TRUE`, draw triangle edges.
#' @param draw.faces If `TRUE`, draw triangles.
#' @param edge.color Edge color specification.
#' @param S Deprecated.
#' @param \dots Additional parameters passed to and from other methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [plot.fm_mesh_2d()]
# @method plot_rgl fm_mesh_2d
#' @export
#' @examples
#' \donttest{
#' if (require(rgl)) {
#'   mesh <- fm_rcdt_2d(globe = 10)
#'   plot_rgl(mesh, col = mesh$loc[, 1])
#' }
#' }
plot_rgl <- function(x, ...) {
  UseMethod("plot_rgl")
}
#' @export
#' @rdname plot_rgl
lines_rgl <- function(x, ..., add = TRUE) {
  UseMethod("lines_rgl")
}
#' @export
#' @rdname plot_rgl
lines_rgl.fm_segm <- function(x, ..., add = TRUE) {
  lines(x, add = add, rgl = TRUE, ...)
}
#' @export
#' @inheritParams plot.fm_mesh_2d
#' @rdname plot_rgl
plot_rgl.fm_mesh_2d <- function(x, col = "white", color.axis = NULL,
                                color.n = 512, color.palette = cm.colors,
                                color.truncate = FALSE, alpha = NULL,
                                lwd = 1, specular = "black",
                                draw.vertices = TRUE,
                                draw.edges = TRUE,
                                draw.faces = TRUE,
                                draw.segments = draw.edges,
                                size = 2,
                                edge.color = rgb(0.3, 0.3, 0.3),
                                t.sub = seq_len(nrow(x$graph$tv)),
                                visibility = "",
                                S = deprecated(),
                                add = FALSE,
                                ...) {
  requireNamespace("rgl", quietly = FALSE)
  mesh <- x

  if (!add) {
    dev <- rgl::open3d()
    rgl::view3d(0, 0, fov = 0)
  } else {
    dev <- rgl::cur3d()
  }

  tv_info <- get_tv_sub(
    tv = mesh$graph$tv,
    loc = mesh$loc,
    t.sub = t.sub,
    visibility = visibility
  )
  tv <- tv_info$tv

  if (draw.vertices) {
    idx <- intersect(unique(as.vector(tv)), mesh$idx$loc)
    rgl::points3d(mesh$loc[idx, , drop = FALSE],
      size = 2 * size, lwd = lwd, color = "blue", ...
    )
  }
  if (draw.segments) {
    if (!is.null(mesh$segm$bnd)) {
      lines_rgl(fm_as_fm(mesh$segm$bnd),
        loc = mesh$loc,
        lwd = lwd + 1,
        ...
      )
    }
    if (!is.null(mesh$segm$int)) {
      lines_rgl(fm_as_fm(mesh$segm$int),
        loc = mesh$loc,
        lwd = lwd + 1,
        ...
      )
    }
  }

  S <- mesh$loc
  TV <- tv

  colors <- (generate_colors(
    col, color.axis, color.n,
    color.palette, color.truncate, alpha
  ))

  tTV <- t(TV)
  Tx <- S[tTV, 1]
  Ty <- S[tTV, 2]
  Tz <- S[tTV, 3]
  if (length(colors$colors) == 1) {
    ## One color
    Tcol <- colors$colors
    Talpha <- colors$alpha
  } else if (length(colors$colors) == nrow(S)) {
    ## One color per vertex
    Tcol <- colors$colors[tTV]
    Talpha <- colors$alpha[tTV]
  } else {
    ## One color per triangle
    stopifnot(length(colors$colors) == nrow(TV))
    Tcol <- colors$colors[t(matrix(rep(seq_len(nrow(TV)), 3), dim(TV)))]
    Talpha <- colors$alpha[t(matrix(rep(seq_len(nrow(TV)), 3), dim(TV)))]
  }
  if (draw.edges) {
    Ev <- rbind(
      cbind(TV[, 1], TV[, 2]),
      cbind(TV[, 2], TV[, 3]),
      cbind(TV[, 3], TV[, 1])
    )
    Ev <- cbind(
      pmin(Ev[, 1], Ev[, 2]),
      pmax(Ev[, 1], Ev[, 2])
    )
    Ev <- unique(Ev)
    if (identical(mesh$manifold, "S2")) {
      # Subdivide
      radius <- mean(rowSums(S^2)^0.5)
      n <- 1 + 2^3
      S2 <- matrix(0.0, nrow(Ev) * n, 3)
      for (k in seq_len(n) - 1L) {
        S2[1L + k + n * (seq_len(nrow(Ev)) - 1L), ] <-
          S[Ev[, 1], ] * (1 - k / (n - 1)) + S[Ev[, 2], ] * k / (n - 1)
      }
      S2 <- S2 / rowSums(S2^2)^0.5 * radius
      Ev2 <- matrix(seq_len(nrow(Ev) * n), nrow(Ev), n, byrow = TRUE)
      # Collect info
      Ev2 <- as.vector(t(cbind(Ev2, NA)))
      Ec <- S2[Ev2, , drop = TRUE]
    } else {
      Ev2 <- as.vector(t(cbind(Ev, NA)))
      Ec <- S[Ev2, , drop = FALSE]
    }
    Ecol <- edge.color
  }
  if (draw.vertices) {
    S_ <- S[unique(as.vector(TV)), , drop = FALSE]
    rgl::points3d(S_, color = "black", ...)
  }
  if (draw.edges) {
    rgl::lines3d(Ec[, 1], Ec[, 2], Ec[, 3], color = Ecol, lwd = lwd, ...)
  }
  if (draw.faces) {
    rgl::triangles3d(Tx, Ty, Tz, color = Tcol, specular = specular, alpha = Talpha, ...)
  }

  return(invisible(dev))
}

## library(geometry)
## S = cbind(x=rnorm(30), y=rnorm(30), z=0)
## TV = delaunayn(S[, 1:2]) # NOTE: inconsistent triangle orders, only for test.
## trimesh(TV, S)
##
## colors = rgb(runif(30), runif(30), runif(30))
## rgl.viewpoint(0, 0, fov=20)
## plot.inla.trimesh(TV, S, colors)

## Ecol = col2rgb(color)/256
## Ecol = Ecol*0.5+(1-0.5)*0 # Rescale towards black
## Ecol = 1-Ecol # Invert
## Ecol = Ecol[, c(2, 3, 1)] # Permute
## Ecol = rgb(Ecol[1,], Ecol[2,], Ecol[3,], maxColorValue = 1)
## Ecol = Ecol[tETV]


#' Draw a triangulation mesh object
#'
#' Plots an [fm_mesh_2d()] object using standard graphics.
#'
#' @param x An [fm_mesh_2d()] object.
#' @param col Color specification.  A single named color, a vector of scalar
#' values, or a matrix of RGB values.  Requires `rgl=TRUE`.
#' @param t.sub Optional triangle index subset to be drawn.
#' @param add If `TRUE`, adds to the current plot instead of starting a
#' new one.
#' @param lwd Line width for triangle edges.
#' @param xlim X-axis limits.
#' @param ylim Y-axis limits.
#' @param main Deprecated.
#' @param size argument `cex` for vertex points.
#' @param draw.vertices If `TRUE`, draw triangle vertices.
#' @param vertex.color Color specification for all vertices.
#' @param draw.edges If `TRUE`, draw triangle edges.
#' @param edge.color Color specification for all edges.
#' @param draw.segments If `TRUE`, draw boundary and interior constraint
#' edges more prominently.
#' @param visibility If "front" only display mesh faces with normal pointing
#' towards the camera.
#' @param asp Aspect ratio for new plots. Default 1.
#' @param \dots Further graphics parameters, interpreted by the respective
#' plotting systems.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [plot.fm_segm()], [plot_rgl.fm_mesh_2d()]
#' @examples
#'
#' mesh <- fm_rcdt_2d(globe = 10)
#' plot(mesh)
#'
#' @rdname plot.fm_mesh_2d
#' @name plot.fm_mesh_2d
#' @export
lines.fm_mesh_2d <- function(x, ..., add = TRUE) {
  plot(x, ..., add = add)
}


get_tv_sub <- function(tv, loc, t.sub, visibility = "front") {
  tv <- tv[t.sub, , drop = FALSE]
  # Filter out away-facing triangles
  if (identical(visibility, "front")) {
    e0 <- loc[tv[, 2], ] - loc[tv[, 1], ]
    e1 <- loc[tv[, 3], ] - loc[tv[, 1], ]
    normal <-
      cbind(
        e0[, 2] * e1[, 3] - e1[, 3] * e1[, 2],
        e0[, 3] * e1[, 1] - e0[, 1] * e1[, 3],
        e0[, 1] * e1[, 2] - e0[, 2] * e1[, 1]
      )
    ok <- normal[, 3] > 0
    tv <- tv[ok, , drop = FALSE]
    t.sub <- t.sub[ok]
  }
  return(list(tv = tv, t.sub = t.sub))
}


#' @rdname plot.fm_mesh_2d
#' @param rgl Deprecated
#' @inheritParams plot.fm_segm
#' @export
#' @examples
#' mesh <- fm_mesh_2d(cbind(0, 1), offset = c(1, 1.5), max.edge = 0.5)
#' plot(mesh)
plot.fm_mesh_2d <- function(
    x,
    col = "white",
    t.sub = seq_len(nrow(x$graph$tv)),
    add = FALSE,
    lwd = 1,
    xlim = range(x$loc[, 1]),
    ylim = range(x$loc[, 2]),
    main = NULL,
    size = 1,
    draw.vertices = FALSE,
    vertex.color = "black",
    draw.edges = TRUE,
    edge.color = rgb(0.3, 0.3, 0.3),
    draw.segments = draw.edges,
    rgl = deprecated(),
    visibility = "front",
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ...) {
  force(t.sub)
  force(xlim)
  force(ylim)
  mesh <- x

  tv_info <- get_tv_sub(
    tv = mesh$graph$tv,
    loc = mesh$loc,
    t.sub = t.sub,
    visibility = visibility
  )
  tv <- tv_info$tv

  idx <- cbind(tv[, c(1:3, 1), drop = FALSE], NA)

  if (draw.edges) {
    TV <- tv
    S <- mesh$loc
    Ev <- rbind(
      cbind(TV[, 1], TV[, 2]),
      cbind(TV[, 2], TV[, 3]),
      cbind(TV[, 3], TV[, 1])
    )
    Ev <- cbind(
      pmin(Ev[, 1], Ev[, 2]),
      pmax(Ev[, 1], Ev[, 2])
    )
    Ev <- unique(Ev)
    if (identical(mesh$manifold, "S2")) {
      # Subdivide
      radius <- mean(rowSums(S^2)^0.5)
      n <- 1 + 2^3
      S2 <- matrix(0.0, nrow(Ev) * n, 3)
      for (k in seq_len(n) - 1L) {
        S2[1L + k + n * (seq_len(nrow(Ev)) - 1L), ] <-
          S[Ev[, 1], ] * (1 - k / (n - 1)) + S[Ev[, 2], ] * k / (n - 1)
      }
      S2 <- S2 / rowSums(S2^2)^0.5 * radius
      Ev2 <- matrix(seq_len(nrow(Ev) * n), nrow(Ev), n, byrow = TRUE)
      # Collect info
      Ev2 <- as.vector(t(cbind(Ev2, NA)))
      Ec <- S2[Ev2, , drop = TRUE]
    } else {
      Ev2 <- as.vector(t(cbind(Ev, NA)))
      Ec <- S[Ev2, , drop = FALSE]
    }
    Ecol <- edge.color
  }


  if (!add) {
    plot(NA, type = "n",
         xlim = xlim, ylim = ylim, asp = asp,
         axes = axes,
         xlab = xlab, ylab = ylab,
         ...)
  }
  if (draw.edges) {
    lines(Ec[, 1], Ec[, 2], type = "l", col = edge.color, lwd = lwd)
  }

  if (draw.vertices) {
    idx <- unique(as.vector(tv))
    points(mesh$loc[idx, , drop = FALSE],
      pch = 20, col = vertex.color, cex = size, ...
    )
    idx <- intersect(idx, mesh$idx$loc)
    points(mesh$loc[idx, , drop = FALSE],
      pch = 20, col = "blue", cex = size, ...
    )
  }
  if (draw.segments) {
    if (!is.null(mesh$segm$bnd)) {
      lines(fm_as_fm(mesh$segm$bnd), mesh$loc, lwd = lwd + 1, ...)
    }
    if (!is.null(mesh$segm$int)) {
      lines(fm_as_fm(mesh$segm$int), mesh$loc, lwd = lwd + 1, ...)
    }
  }
  return(invisible())
}
