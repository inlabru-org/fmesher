# Refine a 2d mesh

Refine an existing mesh

## Usage

``` r
fm_refine(mesh, refine = list(max.edge = 1))
```

## Arguments

- mesh:

  An
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  object

- refine:

  A list of refinement options passed on to
  [fm_rcdt_2d_inla](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)

## Value

A refined `fm_mesh_2d` object

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
fm_dof(fmexample$mesh)
#> [1] 279
fm_dof(fm_refine(fmexample$mesh, refine = list(max.edge = 1)))
#> [1] 332
```
