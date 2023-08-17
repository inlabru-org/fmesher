library(tidyverse)
library(fmesher)
library(INLA)
library(Matrix)
d_ <- c(10, 32, 100, 320, 1000)#, 3200)
prob_ <- c(1, 0.5, 0.2,
           0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.038, 0.036, 0.034, 0.032, 0.03,
           0.02, 0.018, 0.016, 0.014, 0.012,
           0.01, 0.005, 0.002, 0.001)
prob_ <- c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
prob_C_ <- NULL
result <- list()
progressr::handlers("cli")
progressr::with_progress({
  pr <- progressr::progressor(length(d_) * length(prob_))
for (dim in d_) {
  for (prob in prob_) {
    A <- Matrix(rnorm(dim^2), dim, dim)
    B <- (Matrix(runif(dim^2), dim, dim) <= prob) * 1.0
    mean((B%*%t(B)) != 0)
    C <- Diagonal(dim, 1) + (A %*% t(A)) * (B %*% t(B))
    mean(C != 0)
    prob_C_ <- c(prob_C_, mean(C != 0))
    result <- c(result,
                list(bench::mark(fm_qinv(C),
                                 solve(C),
                                 inla.qinv(C),
                                 check = FALSE,
                                 iterations = 1)
                     ))
    pr()
  }
}
})
result <- do.call(rbind, result)
result$prob <- rep(rep(prob_, each = 3), times = length(d_))
result$prob_C <- rep(prob_C_, each = 3)
result$dim <- rep(d_, each = 3 * length(prob_))

ggplot(result) +
  geom_point(aes(dim*prob_C, `itr/sec`, colour = as.character(expression))) +
  scale_x_log10() +
  scale_y_log10()

ggplot(result) +
  geom_point(aes(prob, prob_C, col = factor(dim))) +
  scale_x_log10() +
  scale_y_log10()

