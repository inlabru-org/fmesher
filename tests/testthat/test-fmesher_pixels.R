test_that("fm_pixels sp vs sf", {
  skip_on_cran()
  local_fm_safe_inla()

  mesh <- INLA::inla.mesh.2d(cbind(0,0), offset = 10, max.edge = 1,
                             crs = fm_CRS("longlat_globe"))

  mydata <- sp::SpatialPointsDataFrame(
    mesh$loc,
    data = data.frame(y = rnorm(mesh$n) + 10),
    proj4string = fm_CRS("longlat_globe")
  )

  system.time({
    surface1 <- fm_pixels(mesh, nx = 5, ny = 5, mask = TRUE, format = "sp")
  })

  system.time({
    surface2 <- fm_pixels(mesh, nx = 5, ny = 5, mask = TRUE, format = "sf")
  })

  expect_equal(
    sf::st_coordinates(sf::st_as_sf(surface1)),
    sf::st_coordinates(surface2)
  )

})
