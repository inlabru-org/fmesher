test_that("Discrete integration", {
  domain <- 2:5
  samplers <- 3:7
  ips_ <- data.frame(x = 3:5, weight = rep(1, 3), .block = 1L:3L)

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
  ips_ <- data.frame(
    x = c(3:5, 3.5, 4.5),
    weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
    .block = 1L
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

  # degree = 2
  domain <- fm_mesh_1d(2:5, degree = 2)

  samplers <- c(3, 5)
  ips_ <- data.frame(
    x = c(3:5, 3.5, 4.5),
    weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
    .block = 1L
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
    dplyr::arrange(ips1, .block, time, space),
    dplyr::arrange(ips2[names(ips1)], .block, time, space)
  )
})



# From old ipoints tests

test_that("conversion of SpatialPolygon to integration points when domain is defined via a mesh", {
  ips <- fm_int(fmexample$mesh, samplers = fmexample$boundary_sp[[1]])

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(
    sort(colnames(as.data.frame(ips))),
    sort(c("weight", ".block", "x", "y", "z"))
  )
  expect_equal(sum(ips$weight), 18.33349, tolerance = lowtol)
})

test_that("conversion of whole 2D mesh to integration points", {
  ips <- fm_int(fmexample$mesh, format = "sf")

  expect_s3_class(ips, "sf")
  expect_equal(colnames(ips), c("weight", ".block", "geometry"))
  expect_equal(sum(ips$weight), 64.58135, tolerance = lowtol)

  ips <- fm_int(fmexample$mesh, format = "sp")

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(as.data.frame(ips)), c("weight", ".block", "x", "y", "z"))
  expect_equal(sum(ips$weight), 64.58135, tolerance = lowtol)
})


test_that("Polygon integration with holes", {
  plyA <- sp::SpatialPolygons(list(
    sp::Polygons(
      list(
        sp::Polygon(matrix(c(0, 3, 3, 0, 0, 0, 3, 3) - 1, 4, 2),
          hole = FALSE
        ),
        sp::Polygon(matrix(c(1, 2, 2, 1, 1, 1, 2, 2) - 1, 4, 2),
          hole = TRUE
        )
      ),
      ID = "A"
    )
  ))

  bndA <- fm_as_segm(plyA)
  m <- fmexample$mesh
  ipA1 <- fm_int(m, plyA, int.args = list(
    method = "direct",
    nsub2 = 1
  ))
  ipA2 <- fm_int(m, plyA, int.args = list(
    method = "stable",
    nsub2 = 1
  ))
  ipA3 <- fm_int(m, plyA, int.args = list(method = "direct"))
  ipA4 <- fm_int(m, plyA, int.args = list(method = "stable"))
  ipA1$test <- "A1"
  ipA2$test <- "A2"
  ipA3$test <- "A3"
  ipA4$test <- "A4"

  # if (FALSE) {
  #   require("ggplot2")
  #   pl <- ggplot2::ggplot() +
  #     geom_fm(data = m, alpha = 0) +
  #     inlabru::gg(plyA)
  #   pl
  #
  #   pl +
  #     inlabru::gg(ipA1, mapping = aes(col = weight, size = weight)) +
  #     inlabru::gg(ipA2, mapping = aes(col = weight, size = weight)) +
  #     inlabru::gg(ipA3, mapping = aes(col = weight, size = weight)) +
  #     inlabru::gg(ipA4, mapping = aes(col = weight, size = weight)) +
  #     ggplot2::facet_wrap(vars(test))
  # }

  #   sf::st_area(sf::st_as_sf(plyA))
  # [1] 8.006112

  expect_equal(sum(ipA1$weight), 7.914124, tolerance = lowtol)
  expect_equal(sum(ipA2$weight), 7.914124, tolerance = lowtol)
  expect_equal(sum(ipA3$weight), 8.001529, tolerance = lowtol)
  expect_equal(sum(ipA4$weight), 8.001529, tolerance = lowtol)
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

test_that("flat mesh integration", {
  mesh <- fmexample$mesh

  ips0 <- fm_int(mesh, int.args = list(nsub2 = 0))
  ips9 <- fm_int(mesh, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0$weight), sum(ips9$weight))
})

test_that("sphere and globe mesh integration", {
  skip_on_cran()

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

  suppressWarnings(
    mesh_2 <- fm_rcdt_2d_inla(globe = 1, crs = fm_CRS("globe"))
  )

  ips0_2 <- fm_int(mesh_2, int.args = list(nsub2 = 0))
  ips9_2 <- fm_int(mesh_2, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_2$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_2$weight), 4 * pi * 6370.997^2)

  ips0_3 <- fm_int(mesh_2, int.args = list(nsub2 = 0))
  ips9_3 <- fm_int(mesh_2, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_3$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_3$weight), 4 * pi * 6370.997^2)
})

test_that("flat SpatialPolygons integration", {
  mesh <- fmexample$mesh

  poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(rbind(
    c(-1, -1), c(-1, 1), c(1, 1), c(1, -1)
  ))), ID = "A")))
  poly <- sf::st_as_sf(poly)

  ips0 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 0, method = "direct"))
  ips1 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 1, method = "direct"))
  ips9 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 9, method = "direct"))
  ips19 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 19, method = "direct"))

  # require("ggplot2")
  # ggplot() +
  #   geom_fm(data = mesh) +
  #   geom_sf(aes(size = weight, colour = nsub2), data = cbind(ips0, nsub2 = "0")) +
  #   geom_sf(aes(size = weight, colour = nsub2), data = cbind(ips1, nsub2 = "1")) +
  #   geom_sf(aes(size = weight, colour = nsub2), data = cbind(ips9, nsub2 = "9")) +
  #   geom_sf(aes(size = weight, colour = nsub2), data = cbind(ips19, nsub2 = "19")) +
  #   facet_wrap(~nsub2)

  expect_equal(sum(ips0$weight), 4.055089, tolerance = midtol)
  expect_equal(sum(ips1$weight), 4.02438, tolerance = midtol)
  expect_equal(sum(ips9$weight), 4.00746, tolerance = midtol)
  expect_equal(sum(ips19$weight), 3.999283, tolerance = midtol)
})

test_that("globe polygon integration", {
  skip_on_cran()

  suppressWarnings(
    mesh <- fm_rcdt_2d_inla(globe = 1, crs = fm_CRS("globe"))
  )

  poly <- sp::SpatialPolygons(
    list(sp::Polygons(list(sp::Polygon(rbind(
      c(-45, -45), c(-45, 45), c(45, 45), c(45, -45)
    ))), ID = "A")),
    proj4string = fm_CRS("longlat_globe")
  )

  ips1 <- fm_int(mesh, samplers = poly, int.args = list(nsub = 2))

  expect_equal(
    nrow(ips1),
    9
  )
})
