library(tidyverse)
library(fmesher)
library(INLA)
d <- 50
prob_ <- c(1, 0.5, 0.2,
           0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.038, 0.036, 0.034, 0.032, 0.03,
           0.02, 0.018, 0.016, 0.014, 0.012,
           0.01, 0.005, 0.002, 0.001)
prob_C_ <- NULL
result <- NULL
for (prob in prob_) {
    A <- Matrix(rnorm(d^2), d, d)
  B <- (Matrix(runif(d^2), d, d) <= prob) * 1.0
  mean((B%*%t(B)) != 0)
  C <- Diagonal(d, 1) + (A %*% t(A)) * (B %*% t(B))
  mean(C != 0)
  prob_C_ <- c(prob_C_, mean(C != 0))
  result <- rbind(
    result,
    bench::mark(fm_qinv(C), solve(C), inla.qinv(C), check = FALSE)
  )
}
result$prob <- rep(prob_, each = 3)
result$prob_C <- rep(prob_C_, each = 3)

ggplot(result) +
  geom_point(aes(prob_C, `itr/sec`, colour = as.character(expression))) +
  scale_x_log10() +
  scale_y_log10()

ggplot(result) +
  geom_point(aes(prob, prob_C)) +
  scale_x_log10() +
  scale_y_log10()

