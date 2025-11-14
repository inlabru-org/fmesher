# Recursive curve simplification.

Helper from legacy `INLA::inla.simplify.curve()`

Attempts to simplify a polygonal curve by joining nearly colinear
segments.

Uses a variation of the binary splitting Ramer-Douglas-Peucker
algorithm, with an ellipse of half-width `eps` ellipse instead of a
rectangle, motivated by prediction ellipse for Brownian bridge.

## Usage

``` r
fm_simplify_helper(loc, idx, eps = NULL, eps_rel = NULL)
```

## Arguments

- loc:

  Coordinate matrix.

- idx:

  Index vector into `loc` specifying a polygonal curve.

- eps:

  Absolute straightness tolerance. Default `NULL`, no constraint.

- eps_rel:

  Relative straightness tolerance. Default `NULL`, no constraint.

## Value

An index vector into `loc` specifying the simplified polygonal curve.

## Details

Variation of Ramer-Douglas-Peucker. Uses width epsilon ellipse instead
of rectangle, motivated by prediction ellipse for Brownian bridge.

## See also

Other nonconvex inla legacy support:
[`fm_nonconvex_hull_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull_inla.md),
[`fm_segm_contour_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_segm_contour_helper.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
theta <- seq(0, 2 * pi, length.out = 1000)
loc <- cbind(cos(theta), sin(theta))
idx <- fm_simplify_helper(loc = loc, idx = 1:nrow(loc), eps = 0.01)
print(c(nrow(loc), length(idx)))
#> [1] 1000   33
plot(loc, type = "l")
lines(loc[idx, ], col = "red")
```
