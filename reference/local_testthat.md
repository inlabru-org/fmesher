# Unit test helpers

Local helper functions for package unit tests

## Usage

``` r
local_fm_testthat_assign(x, values, envir = parent.frame())

local_fm_testthat_tolerances(
  tolerances = c(1e-04, 0.01, 0.1),
  envir = parent.frame()
)

local_fm_testthat_setup(envir = parent.frame())
```

## Arguments

- x:

  character; Name of variable to assign to

- values:

  the object to assign to `x`

- envir:

  environment for exit handlers

- tolerances:

  numeric vector of length 3; `[lowtol, midtol, hitol]`

## Value

None

## Functions

- `local_fm_testthat_assign()`: Assign local variable. Useful for easy
  cleanup of global workspace with
  [`withr::deferred_run()`](https://withr.r-lib.org/reference/defer.html)
  when running tests interactively.

- `local_fm_testthat_tolerances()`: Assign test tolerances Assign local
  tolerance variables. Useful for easy cleanup of global workspace with
  [`withr::deferred_run()`](https://withr.r-lib.org/reference/defer.html)
  when running tests interactively.

- `local_fm_testthat_setup()`: Initialise environment for tests. To be
  called either at the top of a testfile, or inside tests.

## Examples

``` r
outer_fun <- function() {
  fun <- function(envir = parent.frame()) {
    local_fm_testthat_assign("local_var_name", 1:4, envir = envir)
  }
  fun()
  local_var_name
}
exists("local_var_name")
#> [1] FALSE
outer_fun()
#> [1] 1 2 3 4
exists("local_var_name")
#> [1] FALSE
```
