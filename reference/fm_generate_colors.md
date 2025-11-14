# Generate text RGB color specifications.

Generates a text RGB color specification matrix based on a color
palette.

## Usage

``` r
fm_generate_colors(
  color,
  color.axis = NULL,
  color.n = 512,
  color.palette = cm.colors,
  color.truncate = FALSE,
  alpha = NULL
)
```

## Arguments

- color:

  `character`, `matrix` or `vector`

- color.axis:

  The min/max limit values for the color mapping.

- color.n:

  The number of colors to use in the color palette.

- color.palette:

  A color palette function.

- color.truncate:

  If `TRUE`, truncate the colors at the color axis limits.

- alpha:

  Transparency/opaqueness values.

## Value

A list with character vector `colors` and numeric vector `alpha`

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
fm_generate_colors(1:4, color.axis = c(1, 4))
#> $colors
#> [1] "#80FFFF" "#D4FFFF" "#FFD4FF" "#FF80FF"
#> 
#> $alpha
#> [1] 1 1 1 1
#> 
```
