# Deprecated functions in fmesher

These functions still attempt to do their job, but will be removed in a
future version.

## Usage

``` r
fm_mesh_components(...)

fm_int_object(...)

fm_sp2segment(...)
```

## Arguments

- ...:

  Usually passed on to other methods

## Functions

- `fm_mesh_components()`: Backwards compatibility for
  [`fm_components()`](https://inlabru-org.github.io/fmesher/reference/fm_components.md),
  deprecated since version `0.4.0.9001`, disabled since `0.6.0`

- `fm_int_object()`: Deprecated function since `0.5.0.9013`; use
  [`new_fm_int()`](https://inlabru-org.github.io/fmesher/reference/new_fm_int.md)
  instead.

- `fm_sp2segment()`: **\[deprecated\]** in favour of
  [`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>
