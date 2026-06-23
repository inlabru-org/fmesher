# Make a collection function space

**\[experimental\]** Collection function spaces. The interface and
object storage model is experimental and may change.

## Usage

``` r
fm_collect(x, ...)
```

## Arguments

- x:

  list of function space objects, such as
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
  all of the same type.

- ...:

  Currently unused

## Value

A `fm_collect` or `fm_collect_list` object. Elements of `fm_collect`:

- fun_spaces:

  `fm_list` of function space objects

- manifold:

  character; manifold type summary, obtained from the function spaces.

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
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md),
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Examples

``` r
m <- fm_collect(list(
  A = fmexample$mesh,
  B = fmexample$mesh
))
m2 <- fm_as_collect(m)
m3 <- fm_as_collect_list(list(m, m))
c(fm_dof(m$fun_spaces[[1]]) + fm_dof(m$fun_spaces[[2]]), fm_dof(m))
#> [1] 584 584
fm_basis(m, loc = tibble::tibble(
  loc = fmexample$loc_sf,
  index = c(1, 1, 2, 2, 1, 2, 2, 1, 1, 2)
), full = TRUE)
#> fm_basis object
#>   Projection matrix (A): 10-by-584
#>   Valid evaluations (ok): 10 out of 10
#>   Additional information: 
fm_basis(m, loc = tibble::tibble(
  loc = rbind(c(0, 0), c(0.1, 0.1)),
  index = c("B", "A")
), full = TRUE)
#> fm_basis object
#>   Projection matrix (A): 2-by-584
#>   Valid evaluations (ok): 2 out of 2
#>   Additional information: 
fm_evaluator(m, loc = tibble::tibble(loc = cbind(0, 0), index = 2))
#> fm_evaluator object
#>   proj:
#>     fm_basis object
#>       Projection matrix (A): 1-by-584
#>       Valid evaluations (ok): 1 out of 1
#>       Additional information: 
#>   Additional evaluator information: 
names(fm_fem(m))
#> [1] "cc" "c0" "c1" "va" "ta" "g1" "g2"
fm_diameter(m)
#> [1] 10.51699
```
