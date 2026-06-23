# Package index

## Mesh construction

Functions for constructing meshes and function spaces.

- [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md)
  : Make a 1D mesh object
- [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  [`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  : Make a 2D mesh object
- [`fm_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_3d.md)
  [`fm_delaunay_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_3d.md)
  : Construct a 3D tetrahedralisation
- [`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)
  [`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)
  [`fm_delaunay_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)
  : Refined Constrained Delaunay Triangulation
- [`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.md)
  **\[experimental\]** : Split triangles of a mesh into subtriangles
- [`fm_subset()`](https://inlabru-org.github.io/fmesher/reference/fm_subset.md)
  **\[experimental\]** : Extract a subset of a mesh
- [`fm_components()`](https://inlabru-org.github.io/fmesher/reference/fm_components.md)
  : Compute connected mesh subsets
- [`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md)
  : Make a lattice object
- [`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md)
  : Lattice grids for N dimensions
- [`fm_hexagon_lattice()`](https://inlabru-org.github.io/fmesher/reference/fm_hexagon_lattice.md)
  **\[experimental\]** : Create hexagon lattice points
- [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  [`fm_segm_join()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  [`fm_segm_split()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  [`fm_is_bnd()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  [`` `fm_is_bnd<-`() ``](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  : Make a spatial segment object
- [`c(`*`<fm_segm>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_segm_list.md)
  [`c(`*`<fm_segm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_segm_list.md)
  [`` `[`( ``*`<fm_segm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_segm_list.md)
  : Methods for fm_segm lists
- [`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
  [`fm_extensions()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
  [`fm_nonconvex_hull_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
  [`fm_nonconvex_hull_sf()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
  : Compute an extension of a spatial object
- [`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md)
  **\[experimental\]** : Recursive curve simplification.
- [`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)
  **\[experimental\]** : Make a tensor product function space
- [`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md)
  **\[experimental\]** : Make a collection function space

## Property methods

Functions for computing or accessing object properties

- [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  : Compute mapping matrix between mesh function space and points

- [`fm_raw_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_raw_basis.md)
  : Basis functions for mesh manifolds

- [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  : Multi-domain integration

- [`new_fm_int()`](https://inlabru-org.github.io/fmesher/reference/new_fm_int.md)
  : Construct integration scheme objects

- [`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  : Compute barycentric coordinates

- [`fm_bary_loc()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_loc.md)
  : Extract Euclidean Sgeometry from Barycentric coordinates

- [`fm_bary_simplex()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_simplex.md)
  : Extract Simplex information for Barycentric coordinates

- [`fm_evaluate()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  [`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  [`fm_evaluator_lattice()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  : Methods for projecting to/from mesh objects

- [`fm_diameter()`](https://inlabru-org.github.io/fmesher/reference/fm_diameter.md)
  : Diameter bound for a geometric object

- [`fm_bbox()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`fm_as_bbox()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`` `[`( ``*`<fm_bbox>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`c(`*`<fm_bbox>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`fm_as_bbox_list()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  : Bounding box class

- [`fm_area()`](https://inlabru-org.github.io/fmesher/reference/fm_area.md)
  : Calculate the area inside segments

- [`fm_sizes()`](https://inlabru-org.github.io/fmesher/reference/fm_sizes.md)
  **\[experimental\]** : fm_sizes

- [`fm_vertices()`](https://inlabru-org.github.io/fmesher/reference/fm_vertices.md)
  :

  Extract vertex locations from an `fm_mesh_2d`

- [`fm_centroids()`](https://inlabru-org.github.io/fmesher/reference/fm_centroids.md)
  :

  Extract triangle centroids from an `fm_mesh_2d`

- [`fm_contains()`](https://inlabru-org.github.io/fmesher/reference/fm_contains.md)
  : Check which mesh triangles are inside a polygon

- [`fm_is_within()`](https://inlabru-org.github.io/fmesher/reference/fm_is_within.md)
  : Query if points are inside a mesh

- [`fm_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  [`fm_manifold_get()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  [`fm_manifold_type()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  [`fm_manifold_dim()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  : Query the mesh manifold type

- [`fm_detect_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_detect_manifold.md)
  [`fm_crs_detect_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_detect_manifold.md)
  : Detect manifold type

- [`fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.md)
  : Function spece degrees of freedom

- [`fm_pixels()`](https://inlabru-org.github.io/fmesher/reference/fm_pixels.md)
  : Generate lattice points covering a mesh

- [`fm_split_lines()`](https://inlabru-org.github.io/fmesher/reference/fm_split_lines.md)
  : Split lines at triangle edges

## Stochastic process methods

Functions for constructing stochastic PDE models

- [`fm_fem()`](https://inlabru-org.github.io/fmesher/reference/fm_fem.md)
  : Compute finite element matrices
- [`fm_matern_precision()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  [`fm_matern_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  [`fm_covariance()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  [`fm_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  **\[experimental\]** : SPDE, GMRF, and Matérn process methods
- [`fm_assess()`](https://inlabru-org.github.io/fmesher/reference/fm_assess.md)
  : Interactive mesh building and diagnostics
- [`fm_qinv()`](https://inlabru-org.github.io/fmesher/reference/fm_qinv.md)
  : Sparse partial inverse

## Matrix helper functions

Functions for handling matrix constructions

- [`fm_block()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  [`fm_block_eval()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  [`fm_block_logsumexp_eval()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  [`fm_block_weights()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  [`fm_block_log_weights()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  [`fm_block_log_shift()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  [`fm_block_prep()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  : Blockwise aggregation matrices
- [`fm_row_kron()`](https://inlabru-org.github.io/fmesher/reference/fm_row_kron.md)
  : Row-wise Kronecker products

## Coordinate system tools

Functions for handling coordinate systems and conversions

- [`` `fm_crs<-`() ``](https://inlabru-org.github.io/fmesher/reference/fm_crs-set.md)
  [`` `fm_crs_oblique<-`() ``](https://inlabru-org.github.io/fmesher/reference/fm_crs-set.md)
  : Assignment operators for crs information objects
- [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  [`fm_crs_oblique()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  [`st_crs(`*`<fm_crs>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  [`` `$`( ``*`<fm_crs>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  [`fm_CRS(`*`<fm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  [`fm_wkt_predef()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  : Obtain coordinate reference system object
- [`fm_crs_is_identical()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_identical.md)
  : Check if two CRS objects are identical
- [`fm_crs_is_null()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_null.md)
  [`is.na(`*`<fm_crs>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_null.md)
  : Check if a crs is NULL or NA
- [`fm_crs_plot()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  [`fm_crs_graticule()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  [`fm_crs_tissot()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  **\[experimental\]** : Plot CRS and fm_crs objects
- [`fm_wkt_is_geocent()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_is_geocent()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_get_ellipsoid_radius()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_get_ellipsoid_radius()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_ellipsoid_radius()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_set_ellipsoid_radius()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`` `fm_ellipsoid_radius<-`() ``](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_set_ellipsoid_radius()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_unit_params()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_get_lengthunit()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_set_lengthunit()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_get_lengthunit()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_set_lengthunit()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_length_unit()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`` `fm_length_unit<-`() ``](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_proj4string()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_tree_projection_type()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_wkt_projection_type()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_projection_type()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  [`fm_crs_bounds()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  : Handling CRS/WKT
- [`fm_detect_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_detect_manifold.md)
  [`fm_crs_detect_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_detect_manifold.md)
  : Detect manifold type
- [`fm_zm()`](https://inlabru-org.github.io/fmesher/reference/fm_zm.md)
  [`fm_zm_input()`](https://inlabru-org.github.io/fmesher/reference/fm_zm.md)
  [`fm_zm_target()`](https://inlabru-org.github.io/fmesher/reference/fm_zm.md)
  **\[experimental\]** : Add or remove Z/M information
- [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.md)
  : Object coordinate transformation
- [`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)
  [`is.na(`*`<fm_CRS>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)
  [`is.na(`*`<inla.CRS>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)
  : Create a coordinate reference system object

## Conversion to and from fmesher classes

Functions for converting non-fmesher objects to fmesher objects.

- [`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md)
  [`fm_as_collect_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md)
  :

  Convert objects to `fm_collect`

- [`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md)
  : Convert objects to fmesher objects

- [`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md)
  [`fm_as_lattice_2d_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md)
  :

  Convert objects to `fm_lattice_2d`

- [`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md)
  [`fm_as_lattice_Nd_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md)
  :

  Convert objects to `fm_lattice_Nd`

- [`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md)
  [`fm_as_mesh_1d_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md)
  :

  Convert objects to `fm_segm`

- [`fm_as_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md)
  [`fm_as_mesh_2d_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md)
  :

  Convert objects to `fm_mesh_2d`

- [`fm_as_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md)
  [`fm_as_mesh_3d_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md)
  :

  Convert objects to `fm_mesh_3d`

- [`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md)
  [`fm_as_segm_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md)
  :

  Convert objects to `fm_segm`

- [`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md)
  : Conversion methods from mesh related objects to sfc

- [`fm_as_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md)
  [`fm_as_tensor_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md)
  :

  Convert objects to `fm_tensor`

- [`fm_bbox()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`fm_as_bbox()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`` `[`( ``*`<fm_bbox>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`c(`*`<fm_bbox>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  [`fm_as_bbox_list()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  : Bounding box class

- [`fm_list()`](https://inlabru-org.github.io/fmesher/reference/fm_list.md)
  [`fm_as_list()`](https://inlabru-org.github.io/fmesher/reference/fm_list.md)
  [`c(`*`<fm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_list.md)
  [`` `[`( ``*`<fm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fm_list.md)
  : Handle lists of fmesher objects

- [`as.triangles3d.fm_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/as.triangles3d.fm_mesh_3d.md)
  : Convert a 3D mesh to a 3D rgl triangulation

## Plotting methods

Functions for plotting.

- [`geom_fm()`](https://inlabru-org.github.io/fmesher/reference/geom_fm.md)
  **\[experimental\]** : ggplot2 geomes for fmesher related objects

- [`lines(`*`<fm_mesh_2d>`*`)`](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.md)
  [`plot(`*`<fm_mesh_2d>`*`)`](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.md)
  : Draw a triangulation mesh object

- [`plot(`*`<fm_segm>`*`)`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.md)
  [`lines(`*`<fm_segm>`*`)`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.md)
  [`plot(`*`<fm_segm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.md)
  [`lines(`*`<fm_segm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.md)
  :

  Draw `fm_segm` objects.

- [`plot_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.md)
  [`lines_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.md)
  : Low level triangulation mesh plotting

- [`fm_crs_plot()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  [`fm_crs_graticule()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  [`fm_crs_tissot()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  **\[experimental\]** : Plot CRS and fm_crs objects

## Printing methods

Functions for printing.

- [`print(`*`<fm_segm>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_segm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_list>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_mesh_2d>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_mesh_3d>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_mesh_1d>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_bbox>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_tensor>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_collect>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_lattice_2d>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_lattice_Nd>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_crs>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  [`print(`*`<fm_CRS>`*`)`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  : Print objects

- [`print(`*`<fm_basis>`*`)`](https://inlabru-org.github.io/fmesher/reference/print.fm_basis.md)
  :

  Print method for `fm_basis`

- [`print(`*`<fm_evaluator>`*`)`](https://inlabru-org.github.io/fmesher/reference/print.fm_evaluator.md)
  :

  Print method for
  [`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)

## Data objects

Data objects

- [`fmexample`](https://inlabru-org.github.io/fmesher/reference/fmexample.md)
  : Example mesh data
- [`fmexample3d`](https://inlabru-org.github.io/fmesher/reference/fmexample3d.md)
  : Example 3D mesh data
- [`fmexample_sp()`](https://inlabru-org.github.io/fmesher/reference/fmexample_sp.md)
  : Add sp data to fmexample

## Direct C++ interface

Direct C++ interface methods, not intended for package users.

- [`fmesher_globe_points()`](https://inlabru-org.github.io/fmesher/reference/fmesher_globe_points.md)
  : Globe points

## Deprecated methods

Deprecated methods
