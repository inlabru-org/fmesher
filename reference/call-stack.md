# Call stack utility functions

Helper functions for displaying call stack information

## Usage

``` r
fm_caller_name(which = 0L, override = NULL)

fm_call_stack(start = 0L, end = 0L, with_numbers = TRUE, ...)

fm_try_callstack(expr)
```

## Arguments

- which:

  The number of frames to go back from the caller

- override:

  character; Overrides the automated function name logic

- start:

  The stack starting point

- end:

  The stack end point

- with_numbers:

  INclude call stack location numbers

- ...:

  Currently unused

- expr:

  An `expression` to evaluate

## Value

`fm_caller_name` returns a string with the the name of a calling
function

`fm_call_stack` returns a character vector

`fm_try_callstack` If successful, returns (invisibly) the value from the
evaluated expression, otherwise an error object with call stack
information attached to the error message.

## Functions

- `fm_call_stack()`:

- `fm_try_callstack()`: Inspired by `berryFunctions::tryStack`

## Examples

``` r
fun <- function() {
  print(fm_caller_name())
  nm <- fm_caller_name()
  print(nm)
}
fun()
#> [1] "print"
#> [1] "fun"
```
