# Defining parameters
lambda <- 1  # Fixed value
lambda1 <- 1  # Fixed value
l <- 4
# Calculate function values


# Create a data frame
kappa_values <- seq(0, l, by = 0.05)
y_values <- sapply(kappa_values,
                   function(kappa) PC_prior_kappa(log_kappa = log(kappa),
                                                  lambda = lambda, lambda1 = lambda1))
df <- data.frame(kappa = kappa_values, y = y_values)

# Plot the prior density of kappa
ggplot(df, aes(x = kappa, y = y)) +
  geom_line() +
  ggtitle(expression(paste("Density of ", kappa)))+
  xlab(expression(kappa)) +
  ylab(expression(paste(pi[kappa], "(", kappa, ")")))+
  theme_minimal()

# Plotting the prior density of v and r
pxl <- expand.grid(v1 = seq(-l, l, length.out = 300), v2 = seq(-l, l, length.out = 300))
pxl$r <-seq(0, 2*l, length.out = 300)
pxl$prior_on_v <- mapply(function(v1, v2) PC_prior_v(v = c(v1,v2), lambda1 = lambda1), pxl$v1, pxl$v2)
pxl$prior_on_r <- mapply(function(r) PC_prior_r(log_r = log(r), lambda1= lambda1), pxl$r)

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

