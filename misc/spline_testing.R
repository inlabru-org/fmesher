library(ggplot2)

z <- sin(seq(-pi/2, pi/2, length.out = 1000))
Z <- cbind(sqrt(1 - z^2), 0, z)
n <- 6
deg <- 3
B1 <- fmesher_spherical_bsplines1(z, n = n, degree = deg, uniform = TRUE)
B2 <- fmesher_spherical_bsplines1(z, n = n, degree = deg, uniform = FALSE)
B3 <- fmesher_spherical_harmonics(Z, max_order = 5, rot_inv = TRUE)
B4 <- fmesher_spherical_harmonics(Z, max_order = 5, rot_int = FALSE)
B3_gsl <- fmesher_spherical_harmonics_gsl(Z, max_order = 5, rot_inv = TRUE)
B4_gsl <- fmesher_spherical_harmonics_gsl(Z, max_order = 5, rot_inv = FALSE)
bench::mark(
  sph_harm_new0 = fmesher_spherical_harmonics(Z, max_order = 5, rot_inv = TRUE),
  sph_harm_gsl0 = fmesher_spherical_harmonics_gsl(Z, max_order = 5, rot_inv = TRUE),
  sph_harm_new = fmesher_spherical_harmonics(Z, max_order = 5, rot_inv = FALSE),
  sph_harm_gsl = fmesher_spherical_harmonics_gsl(Z, max_order = 5, rot_inv = FALSE),
  check = FALSE
)
df <- data.frame(
  z = rep(
    z,
    times = ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B3_gsl) + ncol(B4_gsl)
  ),
  B = c(
    as.vector(B1),
    as.vector(B2),
    as.vector(B3),
    as.vector(B4),
    as.vector(B3_gsl),
    as.vector(B4_gsl)
  ),
  basis = rep(
    c("B-spline", "Spherical", "Spherical (GSL)"),
    length(z) * c(
      ncol(B1) + ncol(B2),
      ncol(B3) + ncol(B4),
      ncol(B3_gsl) + ncol(B4_gsl)
    )
  ),
  type = rep(
    c(
      "Simple",
      "General",
      "Simple",
      "General",
      "Simple",
      "General"
    ),
    length(z) * c(ncol(B1), ncol(B2), ncol(B3), ncol(B4), ncol(B3_gsl), ncol(B4_gsl))
  ),
  index = c(
    rep(seq_len(ncol(B1)), each = nrow(B1)),
    rep(seq_len(ncol(B2)), each = nrow(B2)),
    rep(seq_len(ncol(B3)), each = nrow(B3)),
    rep(seq_len(ncol(B4)), each = nrow(B4)),
    rep(seq_len(ncol(B3_gsl)), each = nrow(B3_gsl)),
    rep(seq_len(ncol(B4_gsl)), each = nrow(B4_gsl))
  )
)
ggplot(df) +
  geom_line(aes((z), B, color = factor(index))) +
  facet_grid(vars(type), vars(basis))

