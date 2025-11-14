# Handle lists of fmesher objects

Methods for constructing and manipulating `fm_list` objects.

## Usage

``` r
fm_list(x, ..., .class_stub = NULL)

fm_as_list(x, ..., .class_stub = NULL)

# S3 method for class 'fm_list'
c(...)

# S3 method for class 'fm_list'
x[i]
```

## Arguments

- x:

  `fm_list` object from which to extract element(s)

- ...:

  Arguments passed to each individual conversion call.

- .class_stub:

  character; class stub name of class to convert each list element to.
  If `NULL`, uses `fm_as_fm` and auto-detects if the resulting list has
  consistent class, and then adds that to the class list. If non-null,
  uses `paste0("fm_as_", .class_stub)` for conversion, and verifies that
  the resulting list has elements consistent with that class.

- i:

  indices specifying elements to extract

## Value

An `fm_list` object, potentially with `fm_{class_stub}_list` added.

## Methods (by generic)

- `c(fm_list)`: The `...` arguments should be coercible to `fm_list`
  objects.

- `[`: Extract sub-list

## Functions

- `fm_list()`: Convert each element of a list, or convert a single
  non-list object and return in a list

- `fm_as_list()`: Convert each element of a list, or convert a single
  non-list object and return in a list

## Examples

``` r
fm_as_list(list(fmexample$mesh, fm_segm_join(fmexample$boundary_fm)))
#> fm_mesh_2d object:
#>   Manifold:  R2
#>   V / E / T: 279 / 805 / 527
#>   Euler char.:   1
#>   Constraints:   Boundary: 29 boundary edges (1 group: 1), Interior: 54 interior edges (1 group: 1)
#>   Bounding box: (-5.331027, 4.061656) x (-3.998161, 5.415609)
#>   Basis d.o.f.:  279
#> 72 boundary edges (1 group: 1)
```
