library(fmesher)
Qinv <- function(Q){
  fmesher:::C_qinv(as(as(as(Q, "Matrix"), "generalMatrix"), "CsparseMatrix"))
}
Qinv_ <- function(Q){
  fmesher:::C_qinv(Q)
}
data <- NULL
progressr::handlers("cli")
progressr::with_progress({
  N_ <- c(1:30, 40, 50, 75, 100)
  N_ <- rev(c(1, 2, 4, 8, 16, 32, 64, 128))
  N_ <- rev(c(1, 2, 4, 8, 16, 32))
  pr <- progressr::progressor(2 * length(N_))
  n <- c()
  for (N in N_) {
    m <- fm_rcdt_2d_inla(globe = N)
    n <- c(n, m$n)
    for (car_order in 1:2) {
      fem <- fm_fem(m)
      if (car_order == 1) {
        Q <- fem$g1 + fem$c0
      } else {
        Q <- fem$g2 + fem$c0
      }
      Q_ <- as(as(Q, "generalMatrix"), "CsparseMatrix")
      timing <-
        bench::mark(S_fm = Qinv(Q),
                    S_fm1 = Qinv_(Q_),
                    S_inla = INLA::inla.qinv(Q),
                    S_inla1 = INLA::inla.qinv(Q, reordering = "amd"),
                    check = FALSE)
      data <- rbind(data,
                    data.frame(n = m$n,
                               N = N,
                               time = as.numeric(timing$median),
                               method = c("fm", "fm1", "inla", "inla1"),
                               car = car_order)
      )
      pr()
    }
  }
})

library(tidyverse)

fit <- list(list(), list())
for (car_order in 1:2) {
  for (method_ in unique(data$method)) {
    fit[[car_order]][[method_]] <- lm(time~I(n/1e3)+I(n^2/1e6)+I(n^3/1e9)+I(n^4/1e12),
                                      data = data %>% filter(method == method_, car==car_order))
    data[data$method == method_ & data$car == car_order, "fit"] <-
      predict(fit[[car_order]][[method_]],
              newdata = data.frame(n = n, N = N_))
  }
}

pl <-
  ggplot(data) +
  geom_point(aes(n, time, col = method)) +
  facet_wrap(vars(car)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_line(aes(n, fit, col = method))
pl
pl <-
  ggplot(data) +
  geom_line(aes(n, time, col = method)) +
  facet_wrap(vars(car)) +
  scale_x_log10() +
  scale_y_log10()
pl

