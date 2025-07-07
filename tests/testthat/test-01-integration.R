test_that("Discrete integration", {
  domain <- 2:5
  samplers <- 3:7
  ips_ <- fm_int_object(
    tibble::tibble(x = 3:5, weight = rep(1, 3), .block = 1L:3L),
    name = "x"
  )

  ips <- fm_int(domain, samplers = samplers)
  expect_identical(ips, ips_)

  domain <- as.character(domain)
  samplers <- as.character(samplers)
  ips <- fm_int(domain, samplers)
  ips_$x <- as.character(ips_$x)
  expect_identical(ips, ips_)

  domain <- factor(domain, levels = domain)
  samplers <- factor(samplers, levels = samplers)
  ips <- fm_int(domain, samplers)
  ips_$x <- factor(ips_$x, levels = domain)
  expect_identical(ips, ips_)
})


test_that("Continuous integration", {
  domain <- fm_mesh_1d(2:5)

  samplers <- c(3, 5)
  ips_ <- fm_int_object(
    tibble::tibble(
      x = c(3:5, 3.5, 4.5),
      weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
      .block = 1L
    ),
    name = "x"
  )
  ips_ <- ips_[order(ips_$x), ]

  ips <- fm_int(domain, samplers = samplers)
  ips <- ips[order(ips$x), ]
  expect_identical(ips, ips_)

  ips_bary_ <- ips_
  ips_bary_$x <- fm_bary(domain, ips_$x)

  ips_bary <- fm_int(domain, samplers = samplers, format = "bary")
  ips_bary_ <- ips_bary_[order(ips_bary_$x$index, ips_bary_$x$where[, 2]), ]
  ips_bary <- ips_bary[order(ips_bary$x$index, ips_bary$x$where[, 2]), ]
  expect_identical(ips_bary, ips_bary_)

  # Check blockwise integration
  samplers <- rbind(c(3, 5), c(2, 4.5))
  ips <- fm_int(domain, samplers = samplers)
  sums <- as.vector(by(ips$weight, ips$.block, sum))
  expect_equal(sums, c(2, 2.5))

  # degree = 2
  domain <- fm_mesh_1d(2:5, degree = 2)

  samplers <- c(3, 5)
  ips_ <- fm_int_object(
    tibble::tibble(
      x = c(3:5, 3.5, 4.5),
      weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
      .block = 1L
    ),
    name = "x"
  )
  ips_ <- ips_[order(ips_$x), ]

  ips <- fm_int(domain, samplers = samplers)
  ips <- ips[order(ips$x), ]
  expect_identical(ips, ips_)

  # Check blockwise integration
  samplers <- rbind(c(3, 5), c(2, 4.5))
  ips <- fm_int(domain, samplers = samplers)
  sums <- as.vector(by(ips$weight, ips$.block, sum))
  expect_equal(sums, c(2, 2.5))
})


test_that("Tensor space integration", {
  mesh_time <- fm_mesh_1d(1:5)
  mesh_space <- fm_mesh_1d(c(0, 5, 10))
  domain <- list(space = mesh_space, time = mesh_time)
  samplers1 <- tibble::tibble(
    time = rbind(c(1, 3), c(2, 4), c(3, 5)),
    space = rbind(c(0, 10), c(0, 5), c(5, 10)),
    weight = c(1, 10, 100)
  )

  ips1 <- fm_int(domain, samplers1)

  samplers2 <- tibble::tibble(
    space = samplers1$space,
    time = samplers1$time,
    weight = samplers1$weight
  )

  ips2 <- fm_int(domain, samplers2)

  expect_equal(sort(names(ips1)), sort(names(ips2)))
  expect_equal(
    dplyr::arrange(
      dplyr::select(ips1, time, space, weight, .block),
      .block, time, space
    ),
    dplyr::arrange(
      dplyr::select(ips2, time, space, weight, .block),
      .block, time, space
    )
  )
})



test_that("Integrating an sf polygon on a mesh domain", {
  ips <- fm_int(fmexample$mesh, samplers = fmexample$boundary_sf[[1]])

  expect_s3_class(ips, "sf")
  expect_equal(
    sort(colnames(ips)),
    sort(c("weight", ".block", ".block_origin", "geometry"))
  )
  expect_equal(sum(ips$weight), 18.339, tolerance = lowtol)

  expect_error(
    fm_int(
      list(geometry = fmexample$mesh),
      samplers = fmexample$boundary_sf[[1]]
    ),
    "Unnamed sampler in the samplers"
  )

  ips <- fm_int(
    list(geometry = fmexample$mesh),
    samplers = list(geometry = fmexample$boundary_sf[[1]])
  )

  expect_s3_class(ips, "sf")
  expect_equal(
    sort(colnames(ips)),
    sort(c("weight", ".block", ".block_origin", "geometry"))
  )
  expect_equal(sum(ips$weight), 18.339, tolerance = lowtol)
})

test_that("Integrating a fm_segm polygon on a mesh domain", {
  ips <- fm_int(fmexample$mesh, samplers = fmexample$boundary_fm[[1]])

  expect_s3_class(ips, "sf")
  expect_equal(
    sort(colnames(ips)),
    sort(c("weight", ".block", ".block_origin", "geometry"))
  )
  expect_equal(sum(ips$weight), 18.339, tolerance = lowtol)
})


test_that("Integrating a SpatialPolygon on a mesh domain", {
  skip_if_not(fm_safe_sp())
  ips <- fm_int(fmexample$mesh, samplers = fmexample_sp()$boundary_sp[[1]])

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(
    sort(colnames(as.data.frame(ips))),
    sort(c("weight", ".block", ".block_origin", "x", "y", "z"))
  )
  expect_equal(sum(ips$weight), 18.339, tolerance = lowtol)
})

test_that("Conversion of whole 2D mesh to integration points", {
  ips1 <- fm_int(fmexample$mesh)
  ips <- fm_int(fmexample$mesh, format = "sf")

  expect_s3_class(ips, "sf")
  expect_equal(
    colnames(ips),
    c("weight", ".block", "geometry", ".block_origin")
  )
  expect_equal(sum(ips$weight), 64.58135, tolerance = lowtol)

  skip_if_not(fm_safe_sp())

  ips <- fm_int(fmexample$mesh, format = "sp")

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(
    names(ips),
    c("weight", ".block", ".block_origin")
  )
  expect_equal(
    sp::coordnames(ips),
    c("x", "y", "z")
  )
  expect_equal(sum(ips$weight), 64.58135, tolerance = lowtol)
})


test_that("Polygon integration with holes", {
  plyA <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(c(0, 3, 3, 0, 0, 0, 0, 3, 3, 0) - 1, 5, 2),
      matrix(c(1, 2, 2, 1, 1, 1, 1, 2, 2, 1) - 1, 5, 2)
    )
  ))

  bndA <- fm_as_segm(plyA)
  m <- fmexample$mesh
  ipA <- fm_int(m,
    plyA,
    int.args = list(method = "direct", nsub2 = 1)
  )

  expect_equal(
    sf::st_area(sf::st_as_sf(plyA)),
    8
  )
  expect_equal(sum(ipA$weight), 7.846134, tolerance = lowtol)
})


test_that("Integration line splitting", {
  mesh <- fmexample$mesh

  segm <- fm_segm(
    loc = rbind(c(-1, 0), c(-1, 1), c(1, 0), c(1, 1)),
    idx = rbind(c(1, 2), c(3, 4)),
    is.bnd = FALSE
  )

  expect_error(
    object = {
      sl <- fm_split_lines(mesh, segm)
    },
    NA
  )

  # Check inlabru issue #63 (problem for single line input), fixed
  segm <- fm_segm(
    loc = rbind(c(-1, 0), c(1, 0)),
    idx = 1:2,
    is.bnd = FALSE
  )
  expect_error(
    object = {
      sl <- fm_split_lines(mesh, segm)
    },
    NA
  )

  # Check if empty input is ok
  segm <- fm_segm(
    loc = NULL,
    idx = integer(0),
    is.bnd = FALSE
  )
  expect_error(
    object = {
      sl <- fm_split_lines(mesh, segm)
    },
    NA
  )
})




# Additional mesh integration tests

test_that("Flat mesh integration", {
  mesh <- fmexample$mesh

  ips0 <- fm_int(mesh, int.args = list(nsub2 = 0))
  ips9 <- fm_int(mesh, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0$weight), sum(ips9$weight))
})

test_that("Sphere and globe mesh integration", {
  mesh <- fm_rcdt_2d_inla(globe = 1)

  ips0 <- fm_int(mesh, int.args = list(nsub2 = 0))
  ips9 <- fm_int(mesh, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0$weight), 4 * pi)
  expect_equal(sum(ips9$weight), 4 * pi)

  mesh_ <- mesh
  mesh_$loc <- mesh_$loc * 1000

  ips0_ <- fm_int(mesh_, int.args = list(nsub2 = 0))
  ips9_ <- fm_int(mesh_, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_$weight), 4 * pi * 1e6)
  expect_equal(sum(ips9_$weight), 4 * pi * 1e6)

  mesh_2 <- fm_rcdt_2d_inla(globe = 1, crs = fm_crs("globe"))

  ips0_2 <- fm_int(mesh_2, int.args = list(nsub2 = 0))
  ips9_2 <- fm_int(mesh_2, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_2$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_2$weight), 4 * pi * 6370.997^2)

  ips0_3 <- fm_int(mesh_2, int.args = list(nsub2 = 0))
  ips9_3 <- fm_int(mesh_2, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_3$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_3$weight), 4 * pi * 6370.997^2)
})


test_that("Globe polygon integration", {
  mesh <- fm_rcdt_2d_inla(globe = 1, crs = fm_crs("globe"))

  poly <- sf::st_sfc(
    sf::st_polygon(
      list(rbind(
        c(-45, -45), c(-45, 45), c(45, 45), c(45, -45), c(-45, -45)
      ))
    ),
    crs = fm_crs("longlat_globe")
  )

  ips1 <- fm_int(mesh, samplers = poly, int.args = list(nsub = 2))

  expect_equal(
    nrow(ips1),
    9
  )
})


test_that("fm_int for linestring coinciding with mesh edges", {
  sc <- 1
  bnd <- fm_segm(
    rbind(
      c(0, 0 * sc),
      c(1, 0 * sc),
      c(1, 1 * sc),
      c(0, 1 * sc)
    ),
    idx = rbind(
      c(1, 2),
      c(2, 3),
      c(3, 4),
      c(4, 1)
    ),
    is.bnd = c(TRUE, TRUE, TRUE, TRUE)
  )
  mesh <- fm_rcdt_2d(boundary = bnd)
  # mesh_sizes <- c(1, 2, 3, 5, 7, 9)
  mesh_sizes <- c(1, 3, 5)
  names(mesh_sizes) <- paste0("mesh", mesh_sizes)
  meshes <- c(
    list(mesh0 = mesh),
    lapply(mesh_sizes, function(n) fm_subdivide(mesh, n = n))
  )

  int_line <- fm_segm(
    rbind(
      c(0.1, 0.1 * sc),
      c(0.95, 0.95 * sc)
    ),
    is.bnd = FALSE
  )

  #  x <- sf::st_linestring(int_line$loc + c(0.01, -0.01, 0.02, -0.03, 0, 0))
  x <- sf::st_linestring(int_line$loc)
  x <- sf::st_as_sf(sf::st_geometry(x))

  int.args <- list(nsub1 = 100)
  ips <- lapply(meshes, fm_int, samplers = x, int.args = int.args)
  (w <- vapply(
    ips,
    function(x) sum(x$weight),
    numeric(1)
  ))
  expect_equal(as.vector(w),
    rep((0.95 - 0.1) * sqrt(1 + sc^2), length(w)),
    tolerance = lowtol
  )

  # ips <- lapply(ips, function(x) x[x$weight > 1e-14, ])
  # (w <- vapply(
  #   ips,
  #   function(x) sum(x$weight),
  #   numeric(1)
  # ))
  # diff(range(w))
  #
  # library(ggplot2)
  # ggplot(
  #   do.call(dplyr::bind_rows,
  #           lapply(c(0, 1, 2, 3, 5, 7, 9, 15),
  #                  function(n) cbind(ips[[paste0("mesh", n)]], n = n)))) +
  #   geom_fm(data = meshes[["mesh3"]]) +
  #   geom_sf(aes(size = weight)) +
  #   scale_size_area() +
  #   geom_sf(data = x, alpha = 0.2, col = "red") +
  #   facet_wrap(~n, nrow = 2)
})

test_that("fm_int for linestring", {
  mesh <- fmexample$mesh
  # mesh_sizes <- c(1, 2, 3, 5, 7, 9)
  mesh_sizes <- c(1, 3)
  names(mesh_sizes) <- paste0("mesh", mesh_sizes)
  meshes <- c(
    list(mesh0 = mesh),
    lapply(mesh_sizes, function(n) fm_subdivide(mesh, n = n))
  )

  x <- sf::st_linestring(fm_as_segm(fmexample$boundary_sf[[1]])$loc)
  x <- sf::st_as_sf(sf::st_geometry(x))

  ips <- lapply(meshes, fm_int, samplers = x)
  (w <- vapply(
    ips,
    function(x) sum(x$weight),
    numeric(1)
  ))
  expect_equal(as.vector(w),
    rep(16.3259194526, length(w)),
    tolerance = lowtol
  )

  # ips <- lapply(ips, function(x) x[x$weight > 1e-14, ])
  # (w <- vapply(
  #   ips,
  #   function(x) sum(x$weight),
  #   numeric(1)
  # ))
  # diff(range(w))
  #
  # library(ggplot2)
  # ggplot(
  #   do.call(dplyr::bind_rows,
  #           lapply(c(0, 1, 2, 3, 5, 7, 9, 15),
  #                  function(n) cbind(ips[[paste0("mesh", n)]], n = n)))) +
  #   geom_fm(data = meshes[["mesh3"]]) +
  #   geom_sf(aes(size = weight)) +
  #   geom_sf(data = fmexample$boundary_sf[[1]], alpha = 0.2) +
  #   scale_size_area() +
  #   facet_wrap(~n, nrow = 2) +
  #   xlim(c(-1.8, -1.4)) +
  #   ylim(c(-1.6, -1.2))
})


test_that("Block integration has correct result order", {
  # Single domain
  samplers <- sf::st_sf(
    geometry = rep(fmexample$boundary_sf[[1]], 10),
    weight = 1:10
  )

  ips0 <- fm_int(fmexample$mesh,
    samplers = fmexample$boundary_sf[[1]],
    int.args = list(nsub2 = 1)
  )
  ips <- fm_int(fmexample$mesh,
    samplers = samplers,
    int.args = list(nsub2 = 1)
  )

  expect_s3_class(ips, "sf")

  vals <- fm_block_eval(
    block = ips$.block,
    weights = ips$weight,
    n_block = 10,
    rescale = FALSE,
    values = rep(1, nrow(ips))
  )

  expect_equal(length(vals), nrow(samplers))
  expect_equal(vals, (1:10) * sum(ips0$weight))

  # Multi-domain, at least one with >= 10 blocks
  samplers <- sf::st_sf(
    geometry = rep(fmexample$boundary_sf[[1]], 10),
    weight = 1:10
  )

  ips <- fm_int(
    list(
      geometry = fmexample$mesh,
      time = 1:2
    ),
    samplers = list(
      samplers,
      data.frame(time = 1:2, weight = c(1, 100))
    ),
    int.args = list(nsub2 = 1)
  )

  expect_s3_class(ips, "sf")

  vals <- fm_block_eval(
    block = ips$.block,
    weights = ips$weight,
    n_block = 20,
    rescale = FALSE,
    values = rep(1, nrow(ips))
  )

  expect_equal(length(vals), nrow(samplers) * 2)
  expect_equal(vals, c((1:10), (1:10) * 100) * sum(ips0$weight))
})
