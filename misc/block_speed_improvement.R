library(fmesher)
fm_block_log_shift_old <- function(block = NULL,
                               log_weights = NULL,
                               n_block = NULL) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      force_log = TRUE
    )
  if (length(info$block) == 0) {
    return(0.0)
  }
  block_k <- sort(unique(info$block))
  shift <- numeric(info$n_block)
  if (!is.null(info$log_weights)) {
    shift[block_k] <-
      vapply(
        block_k,
        function(k) {
          max(info$log_weights[info$block == k])
        },
        0.0
      )
  }

  shift
}

fm_block_log_shift_new <- function(block = NULL,
                                   log_weights = NULL,
                                   n_block = NULL) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      force_log = TRUE
    )
  if (length(info$block) == 0) {
    return(0.0)
  }
  block_k <- sort(unique(info$block))
  shift <- numeric(info$n_block)
  if (!is.null(info$log_weights)) {
    shift[block_k] <-
      stats::aggregate(info$log_weights,
                       by = list(block=info$block),
                       FUN = max, simplify = TRUE)$x
  }

  shift
}


n_block <- 800
block <- rep(seq_len(n_block), times = 40)
log_weights <- runif(length(block))

print(bench::mark(old = fm_block_log_shift_old(block = block, log_weights = log_weights,
                                       n_block = n_block),
            new = fm_block_log_shift_new(block = block, log_weights = log_weights,
                                       n_block = n_block)))
n_block <- 16000
block <- rep(seq_len(n_block), 1)
log_weights <- runif(length(block))

print(bench::mark(old = fm_block_log_shift_old(block = block, log_weights = log_weights,
                                       n_block = n_block),
            new = fm_block_log_shift_new(block = block, log_weights = log_weights,
                                       n_block = n_block)))
n_block <- 2
block <- rep(seq_len(n_block), each = 8000)
log_weights <- runif(length(block))

print(bench::mark(old = fm_block_log_shift_old(block = block, log_weights = log_weights,
                                       n_block = n_block),
            new = fm_block_log_shift_new(block = block, log_weights = log_weights,
                                       n_block = n_block)))
# profvis::profvis(for(loop in seq_len(1000)){
#   fm_block_log_shift_new(block = block, log_weights = log_weights,
#                          n_block = n_block)
# })
library(tidyverse)
times <- NULL
for(n_block in c(1,2,4,8,16,32,64,128,256,512,1024)){
  for(N in 1024*c(1,4,16,64)){
    print(c(n_block, N))
    block <- rep(seq_len(n_block), each = N/n_block)
    log_weights <- runif(length(block))
    times <- bind_rows(
      times,
      bind_cols(
        n_block = n_block,
        N=N,
        bench::mark(
          old = fm_block_log_shift_old(block = block, log_weights = log_weights,
                                       n_block = n_block),
          new = fm_block_log_shift_new(block = block, log_weights = log_weights,
                                       n_block = n_block),
          pkg = fm_block_log_shift(block = block, log_weights = log_weights,
                                       n_block = n_block),
          filter_gc = FALSE)))
  }
}
times$method <- as.character(times$expression)
ggplot(times) + geom_line(aes(x=n_block, y=as.numeric(median), col=method)) +
  facet_wrap(~ N) + scale_y_log10()

