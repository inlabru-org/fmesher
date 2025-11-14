# Make a tensor product function space

**\[experimental\]** Tensor product function spaces. The interface and
object storage model is experimental and may change.

## Usage

``` r
fm_tensor(x, ...)
```

## Arguments

- x:

  list of function space objects, such as
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md).

- ...:

  Currently unused

## Value

A `fm_tensor` or `fm_tensor_list` object. Elements of `fm_tensor`:

- fun_spaces:

  `fm_list` of function space objects

- manifold:

  character; manifold type summary. Regular subset of Rd "Rd", if all
  function spaces have type "R", torus connected "Td" if all function
  spaces have type "S", and otherwise "Md" In all cases, `d` is the sum
  of the manifold dimensions of the function spaces.

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
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md)

## Examples

``` r
m <- fm_tensor(list(
  space = fmexample$mesh,
  time = fm_mesh_1d(1:5)
))
m2 <- fm_as_tensor(m)
m3 <- fm_as_tensor_list(list(m, m))
c(fm_dof(m$fun_spaces$space) * fm_dof(m$fun_spaces$time), fm_dof(m))
#> [1] 1395 1395
str(fm_evaluator(m, loc = list(space = cbind(0, 0), time = 2.5)))
#> List of 1
#>  $ proj:List of 2
#>   ..$ A :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:6] 0 0 0 0 0 0
#>   .. .. ..@ p       : int [1:1396] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..@ Dim     : int [1:2] 1 1395
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : NULL
#>   .. .. ..@ x       : num [1:6] 0.2394 0.2112 0.0494 0.2394 0.2112 ...
#>   .. .. ..@ factors : list()
#>   ..$ ok: logi TRUE
#>   ..- attr(*, "class")= chr "fm_basis"
#>  - attr(*, "class")= chr "fm_evaluator"
str(fm_basis(m, loc = list(space = cbind(0, 0), time = 2.5)))
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:6] 0 0 0 0 0 0
#>   ..@ p       : int [1:1396] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..@ Dim     : int [1:2] 1 1395
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:6] 0.2394 0.2112 0.0494 0.2394 0.2112 ...
#>   ..@ factors : list()
str(fm_fem(m))
#> List of 3
#>  $ cc:Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:1395] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ j       : int [1:1395] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ Dim     : int [1:2] 1395 1395
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   .. ..@ x       : num [1:1395] 0.171 0.271 0.213 0.282 0.22 ...
#>   .. ..@ factors : list()
#>  $ g1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:11677] 0 1 25 257 279 0 1 2 226 257 ...
#>   .. ..@ p       : int [1:1396] 0 5 11 16 22 28 33 39 45 51 ...
#>   .. ..@ Dim     : int [1:2] 1395 1395
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   .. ..@ x       : num [1:11677] 1.337 -0.174 -0.113 -0.708 -0.343 ...
#>   .. ..@ factors : list()
#>  $ g2:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:42761] 0 1 2 24 25 225 226 229 233 245 ...
#>   .. ..@ p       : int [1:1396] 0 17 36 53 73 91 107 127 147 167 ...
#>   .. ..@ Dim     : int [1:2] 1395 1395
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   .. ..@ x       : num [1:42761] 11.9277 -1.7904 0.2465 0.0254 -1.1944 ...
#>   .. ..@ factors : list()
```
