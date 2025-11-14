# Internal WKT handling

Conversion between WKT and a tree representation

## Usage

``` r
fm_wkt_as_wkt_tree(x, ...)

fm_wkt_tree_as_wkt(x, pretty = FALSE, ...)

fm_wkt_tree_get_item(x, item, duplicate = 1)

fm_wkt_tree_set_item(x, item_tree, duplicate = 1)
```

## Arguments

- x:

  A WKT2 string, or a `wkt_tree` list structure

- ...:

  Unused

- pretty:

  logical; If TRUE, use pretty formatting. Default: FALSE

- item:

  character vector with item labels identifying a parameter item entry.

- duplicate:

  For items that have more than one match, `duplicate` indicates the
  index number of the desired version. Default: 1

- item_tree:

  An item tree identifying a parameter item entry

## Value

A hierarchical list, describing WKT information as a tree

`fm_wkt_tree_as_wkt` character; the WKT corresponding to the tree.

`fm_wkt_tree_get_item` returns the value of an item found in the tree

`fm_wkt_tree_set_item` returns the modified tree

## Examples

``` r
str(fm_wkt_as_wkt_tree(fm_crs("longlat_norm")$wkt))
#> List of 2
#>  $ label : chr "GEOGCRS"
#>  $ params:List of 6
#>   ..$ : chr "\"unknown\""
#>   ..$ :List of 2
#>   .. ..$ label : chr "DATUM"
#>   .. ..$ params:List of 2
#>   .. .. ..$ : chr "\"unknown\""
#>   .. .. ..$ :List of 2
#>   .. .. .. ..$ label : chr "ELLIPSOID"
#>   .. .. .. ..$ params:List of 4
#>   .. .. .. .. ..$ : chr "\"unknown\""
#>   .. .. .. .. ..$ : chr "1"
#>   .. .. .. .. ..$ : chr "0"
#>   .. .. .. .. ..$ :List of 2
#>   .. .. .. .. .. ..$ label : chr "LENGTHUNIT"
#>   .. .. .. .. .. ..$ params:List of 3
#>   .. .. .. .. .. .. ..$ : chr "\"metre\""
#>   .. .. .. .. .. .. ..$ : chr "1"
#>   .. .. .. .. .. .. ..$ :List of 2
#>   .. .. .. .. .. .. .. ..$ label : chr "ID"
#>   .. .. .. .. .. .. .. ..$ params:List of 2
#>   .. .. .. .. .. .. .. .. ..$ : chr "\"EPSG\""
#>   .. .. .. .. .. .. .. .. ..$ : chr "9001"
#>   ..$ :List of 2
#>   .. ..$ label : chr "PRIMEM"
#>   .. ..$ params:List of 3
#>   .. .. ..$ : chr "\"Reference meridian\""
#>   .. .. ..$ : chr "0"
#>   .. .. ..$ :List of 2
#>   .. .. .. ..$ label : chr "ANGLEUNIT"
#>   .. .. .. ..$ params:List of 3
#>   .. .. .. .. ..$ : chr "\"degree\""
#>   .. .. .. .. ..$ : chr "0.0174532925199433"
#>   .. .. .. .. ..$ :List of 2
#>   .. .. .. .. .. ..$ label : chr "ID"
#>   .. .. .. .. .. ..$ params:List of 2
#>   .. .. .. .. .. .. ..$ : chr "\"EPSG\""
#>   .. .. .. .. .. .. ..$ : chr "9122"
#>   ..$ :List of 2
#>   .. ..$ label : chr "CS"
#>   .. ..$ params:List of 2
#>   .. .. ..$ : chr "ellipsoidal"
#>   .. .. ..$ : chr "2"
#>   ..$ :List of 2
#>   .. ..$ label : chr "AXIS"
#>   .. ..$ params:List of 4
#>   .. .. ..$ : chr "\"longitude\""
#>   .. .. ..$ : chr "east"
#>   .. .. ..$ :List of 2
#>   .. .. .. ..$ label : chr "ORDER"
#>   .. .. .. ..$ params:List of 1
#>   .. .. .. .. ..$ : chr "1"
#>   .. .. ..$ :List of 2
#>   .. .. .. ..$ label : chr "ANGLEUNIT"
#>   .. .. .. ..$ params:List of 3
#>   .. .. .. .. ..$ : chr "\"degree\""
#>   .. .. .. .. ..$ : chr "0.0174532925199433"
#>   .. .. .. .. ..$ :List of 2
#>   .. .. .. .. .. ..$ label : chr "ID"
#>   .. .. .. .. .. ..$ params:List of 2
#>   .. .. .. .. .. .. ..$ : chr "\"EPSG\""
#>   .. .. .. .. .. .. ..$ : chr "9122"
#>   ..$ :List of 2
#>   .. ..$ label : chr "AXIS"
#>   .. ..$ params:List of 4
#>   .. .. ..$ : chr "\"latitude\""
#>   .. .. ..$ : chr "north"
#>   .. .. ..$ :List of 2
#>   .. .. .. ..$ label : chr "ORDER"
#>   .. .. .. ..$ params:List of 1
#>   .. .. .. .. ..$ : chr "2"
#>   .. .. ..$ :List of 2
#>   .. .. .. ..$ label : chr "ANGLEUNIT"
#>   .. .. .. ..$ params:List of 3
#>   .. .. .. .. ..$ : chr "\"degree\""
#>   .. .. .. .. ..$ : chr "0.0174532925199433"
#>   .. .. .. .. ..$ :List of 2
#>   .. .. .. .. .. ..$ label : chr "ID"
#>   .. .. .. .. .. ..$ params:List of 2
#>   .. .. .. .. .. .. ..$ : chr "\"EPSG\""
#>   .. .. .. .. .. .. ..$ : chr "9122"
```
