Prior_kappa <- function(kappa, lambda, lambda1) {
  c <- 2 * sqrt(3*pi)
  fact <- exp(-((kappa * lambda) / c)) * lambda * lambda1 / (kappa * lambda + lambda1)^2
  return(fact*(1+(kappa*lambda+lambda1)/c))
}
f_r <- function(r){
  return(sqrt((1+3 * cosh(2 * r))/48/pi))
}
f0_r <- f_r(0)
Prior_f <- function(ff, lambda1) {
  if (ff >= f0_r) {
    return(lambda1 * exp(-lambda1 * (ff - f0_r)))
  } else {
    return(0)
  }
}

Df_r <- function(r) {
  return((sqrt(3 / pi) * sinh(2 * r)) / (4 * sqrt(1 + 3 * cosh(2 * r))))
}

Prior_r <- function(r, lambda1) {
  if (r >= 0) {
    ff <- f_r(r)
    return(Prior_f(ff, lambda1) * Df_r(r))
  } else {
    return(0)
  }

}
Prior_v <- function(v1, v2, lambda1) {
  norm_v <- sqrt(v1^2+v2^2)  # Calculate the norm of vector v
  return(Prior_r(norm_v, lambda1) / (2 * pi * norm_v))
}
# Generate data
kappa_values <- seq(0, 10, by = 0.1)
lambda <- 1  # Fixed value
lambda1 <- 1  # Fixed value

# Calculate function values
y_values <- sapply(kappa_values, function(kappa) Prior_kappa(kappa, lambda, lambda1))

# Create a data frame
df <- data.frame(kappa = kappa_values, y = y_values)

# Plot the prior density of kappa
ggplot(df, aes(x = kappa, y = y)) +
  geom_line() +
  ggtitle(expression(paste("Density of ", kappa)))+
  xlab(expression(kappa)) +
  ylab(expression(paste(pi[kappa], "(", kappa, ")")))+
  theme_minimal()

# Plotting the prior denisty of v and r
l <- 4
lambda1 <- 1
pxl <- expand.grid(v1 = seq(-l, l, length.out = 300), v2 = seq(-l, l, length.out = 300))
pxl$r <-seq(0, 2*l, length.out = 300)
pxl$prior_on_v <- mapply(function(v1, v2) Prior_v(v1, v2, lambda1), pxl$v1, pxl$v2)
pxl$prior_on_r <- mapply(function(v1) Prior_r(v1, lambda1), pxl$r)

# Plot
library(ggplot2)
ggplot(pxl, aes(x = v1, y = v2, fill = prior_on_v)) +
  geom_tile() +
  ggtitle(expression(paste("Prior on ", v))) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  labs(x = "v1 Coordinate", y = "v2 Coordinate")

# Assuming pxl is a data frame with columns 'prior_on_r' and 'v1'
ggplot(pxl, aes(x = r, y = prior_on_r)) +
  geom_line() +
  ggtitle(expression(paste("Prior on ", r,"=|v|"))) +
  xlab("r") +
  ylab(expression(paste(pi[r], "(", r, ")")))+
  theme_minimal()

