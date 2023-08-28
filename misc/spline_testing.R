library(ggplot2)

z <- seq(-1, 1, length.out = 1000)
n <- 6
deg <- 3
B1 <- fmesher_spherical_bsplines1(z, n = n, degree = deg, uniform = TRUE)
B2 <- fmesher_spherical_bsplines1(z, n = n, degree = deg, uniform = FALSE)
ggplot(data.frame(
  z = rep(z, times = ncol(B1)),
  B1 = as.vector(B1),
  B2 = as.vector(B2),
  basis = factor(rep(seq_len(ncol(
    B1
  )), each = nrow(B1)))
)) +
  geom_line(aes(asin(z), B1, color = basis, lty = "uniform")) +
  geom_line(aes(asin(z), B2, color = basis, lty = "non-uniform"))
