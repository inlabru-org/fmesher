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
block <- rep(seq_len(n_block), 40)
log_weights <- runif(length(block))

print(bench::mark(old = fm_block_log_shift_old(block = block, log_weights = log_weights,
                                       n_block = n_block),
            new = fm_block_log_shift_new(block = block, log_weights = log_weights,
                                       n_block = n_block)))
