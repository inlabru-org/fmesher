# Recursive curve simplification.

**\[experimental\]** Simplifies polygonal curve segments by joining
nearly co-linear segments.

Uses a variation of the binary splitting Ramer-Douglas-Peucker
algorithm, with an ellipse of half-width `eps` ellipse instead of a
rectangle, motivated by prediction ellipse for Brownian bridge.

## Usage

``` r
fm_simplify(x, eps = NULL, eps_rel = NULL, ...)
```

## Arguments

- x:

  An
  [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  object.

- eps:

  Absolute straightness tolerance. Default `NULL`, no constraint.

- eps_rel:

  Relative straightness tolerance. Default `NULL`, no constraint.

- ...:

  Currently unused.

## Value

The simplified
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
object.

## Details

Variation of Ramer-Douglas-Peucker. Uses width epsilon ellipse instead
of rectangle, motivated by prediction ellipse for Brownian bridge.

## References

Ramer, Urs (1972). "An iterative procedure for the polygonal
approximation of plane curves". *Computer Graphics and Image
Processing*. **1** (3): 244–256.
[doi:10.1016/S0146-664X(72)80017-0](https://doi.org/10.1016/S0146-664X%2872%2980017-0)

Douglas, David; Peucker, Thomas (1973). "Algorithms for the reduction of
the number of points required to represent a digitized line or its
caricature". *The Canadian Cartographer*. **10** (2): 112–122.
[doi:10.3138/FM57-6770-U75U-7727](https://doi.org/10.3138/FM57-6770-U75U-7727)

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
[`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md),
[`fm_as_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md),
[`fm_as_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md),
[`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md),
[`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md),
[`fm_as_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md),
[`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md),
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md),
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
theta <- seq(0, 2 * pi, length.out = 1000)
(segm <- fm_segm(cbind(cos(theta), sin(theta)),
  idx = seq_along(theta)
))
#> fm_segm object:
#>   999 boundary edges
#>   Bounding box = (-0.9999951, 1.0000000) x (-0.9999988, 0.9999988) x (0,0)
(segm1 <- fm_simplify(segm, eps_rel = 0.1))
#> fm_segm object:
#>   16 boundary edges (1 group: 0)
#>   Bounding box = (-0.9999951, 1.0000000) x (-0.9999988, 0.9999889) x (0,0)
(segm2 <- fm_simplify(segm, eps_rel = 0.2))
#> fm_segm object:
#>   8 boundary edges (1 group: 0)
#>   Bounding box = (-0.9999951, 1.0000000) x (-0.9999988, 0.9999889) x (0,0)
plot(segm)
lines(segm1, col = 2)
lines(segm2, col = 3)


(segm <- fm_segm(cbind(theta, sin(theta * 4)),
  idx = seq_along(theta)
))
#> fm_segm object:
#>   999 boundary edges
#>   Bounding box = (0.000000,6.283185) x (-0.9999988, 0.9999988) x (0,0)
(segm1 <- fm_simplify(segm, eps_rel = 0.1))
#> fm_segm object:
#>   73 boundary edges (1 group: 0)
#>   Bounding box = (0.000000,6.283185) x (-0.9999988, 0.9999988) x (0,0)
(segm2 <- fm_simplify(segm, eps_rel = 0.2))
#> fm_segm object:
#>   17 boundary edges (1 group: 0)
#>   Bounding box = (0.000000,6.283185) x (-0.9995538, 0.9994549) x (0,0)
plot(segm)
lines(segm1, col = 2)
lines(segm2, col = 3)
```
