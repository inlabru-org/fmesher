# Blockwise aggregation matrices

Creates an aggregation matrix for blockwise aggregation, with optional
weighting.

## Usage

``` r
fm_block(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL
)

fm_block_eval(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL,
  values = NULL
)

fm_block_logsumexp_eval(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL,
  values = NULL,
  log = TRUE
)

fm_block_weights(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL
)

fm_block_log_weights(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL
)

fm_block_log_shift(block = NULL, log_weights = NULL, n_block = NULL)

fm_block_prep(
  block = NULL,
  log_weights = NULL,
  weights = NULL,
  n_block = NULL,
  values = NULL,
  n_values = NULL,
  force_log = FALSE
)
```

## Arguments

- block:

  integer vector; block information. If `NULL`, `rep(1L, block_len)` is
  used, where `block_len` is determined by `length(log_weights)))` or
  `length(weights)))`. A single scalar is also repeated to a vector of
  corresponding length to the weights.

  Note: from version `0.2.0.9017` to `0.4.0.9005`, 'character' input was
  converted to integer with `as.integer(factor(block))`. As this could
  lead to unintended ordering of the output, this is no longer allowed.

- weights:

  Optional weight vector

- log_weights:

  Optional `log(weights)` vector. Overrides `weights` when non-NULL.

- rescale:

  logical; If `TRUE`, normalise the weights by `sum(weights)` or
  `sum(exp(log_weights))` within each block. Default: `FALSE`

- n_block:

  integer; The number of conceptual blocks. Only needs to be specified
  if it's larger than `max(block)`, or to keep the output of consistent
  size for different inputs.

- values:

  Vector to be blockwise aggregated

- log:

  If `TRUE` (default), return log-sum-exp. If `FALSE`, return sum-exp.

- n_values:

  When supplied, used instead of `length(values)` to determine the value
  vector input length.

- force_log:

  When `FALSE` (default), passes either `weights` and `log_weights` on,
  if provided, with `log_weights` taking precedence. If `TRUE`, forces
  the computation of `log_weights`, whether given in the input or not.

## Value

A (sparse) matrix

## Functions

- `fm_block()`: A (sparse) matrix of size `n_block` times
  `length(block)`.

- `fm_block_eval()`: Evaluate aggregation. More efficient alternative to
  to `as.vector(fm_block(...) %*% values)`.

- `fm_block_logsumexp_eval()`: Evaluate log-sum-exp aggregation. More
  efficient and numerically stable alternative to to
  `log(as.vector(fm_block(...) %*% exp(values)))`.

- `fm_block_weights()`: Computes (optionally) blockwise renormalised
  weights

- `fm_block_log_weights()`: Computes (optionally) blockwise renormalised
  log-weights

- `fm_block_log_shift()`: Computes shifts for stable blocked
  log-sum-exp. To compute \\\log(\sum\_{i; \textrm{block}\_i=k}
  \exp(v_i) w_i)\\ for each block `k`, first compute combined values and
  weights, and a shift:

      w_values <- values + fm_block_log_weights(block, log_weights = log_weights)
      shift <- fm_block_log_shift(block, log_weights = w_values)

  Then aggregate the values within each block:

      agg <- aggregate(exp(w_values - shift[block]),
                       by = list(block = block),
                       \(x) log(sum(x)))
      agg$x <- agg$x + shift[agg$block]

  The implementation uses a faster method:

      as.vector(
        Matrix::sparseMatrix(
          i = block,
          j = rep(1L, length(block)),
          x = exp(w_values - shift[block]),
          dims = c(n_block, 1))
      ) + shift

- `fm_block_prep()`: Helper function for preparing `block`, `weights`,
  and `log_weights`, `n_block` inputs.

## Examples

``` r
block <- rep(1:2, 3:2)
fm_block(block)
#> 2 x 5 sparse Matrix of class "dgCMatrix"
#>               
#> [1,] 1 1 1 . .
#> [2,] . . . 1 1
fm_block(block, rescale = TRUE)
#> 2 x 5 sparse Matrix of class "dgCMatrix"
#>                                           
#> [1,] 0.3333333 0.3333333 0.3333333 .   .  
#> [2,] .         .         .         0.5 0.5
fm_block(block, log_weights = -2:2, rescale = TRUE)
#> 2 x 5 sparse Matrix of class "dgCMatrix"
#>                                                       
#> [1,] 0.09003057 0.2447285 0.665241 .         .        
#> [2,] .          .         .        0.2689414 0.7310586
fm_block_eval(
  block,
  weights = 1:5,
  rescale = TRUE,
  values = 11:15
)
#> [1] 12.33333 14.55556
fm_block_logsumexp_eval(
  block,
  weights = 1:5,
  rescale = TRUE,
  values = log(11:15),
  log = FALSE
)
#> [1] 12.33333 14.55556
```
