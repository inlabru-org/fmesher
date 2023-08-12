#ifdef FMESHER_WITH_R

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Order must be RcppFmesher, Rcpp to allow the Rcpp classes
// to find the fmesher types
#include "RcppFmesher.h"
#include "Rcpp.h"

#include "fmesher.h"
#include "fmesher_helpers.h"

using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;

using fmesh::constrListT;
using fmesh::constrMetaT;
using fmesh::constrT;
using fmesh::Dart;
using fmesh::DartList;
using fmesh::DartPair;
using fmesh::Int3;
using fmesh::Int3Raw;
using fmesh::IOHelper;
using fmesh::IOHelperM;
using fmesh::IOHelperSM;
using fmesh::Matrix;
using fmesh::Matrix3double;
using fmesh::Matrix3int;
using fmesh::Matrix1int;
using fmesh::MatrixC;
using fmesh::Mesh;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::SparseMatrix;
using fmesh::TriangleLocator;
using fmesh::Vector3;
using fmesh::vertexListT;

#ifdef FMESHER_WITH_EIGEN
template <class T> using EigenM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <class T> using EigenM1 = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <class T> using EigenSM = Eigen::SparseMatrix<T>;
template <class T> using EigenMM = Eigen::Map<EigenM<T>>;
template <class T> using EigenMM1 = Eigen::Map<EigenM1<T>>;
template <class T> using EigenMSM = Eigen::Map<EigenSM<T>>;
#endif

// const bool useVT = true;
// const bool useTTi = true;

template <typename T>
bool Rcpp_is_element(const Rcpp::List& list, std::string name) {
  if (!list.containsElementNamed(name.c_str()))
    return false;

  if (Rf_isNull(list[name.c_str()]))
    return false;

  return Rcpp::is<T>(list[name.c_str()]);
}


class Options {
public:
  double cutoff;
  double sphere_tolerance;
  int cet_sides;
  double cet_margin;
  double rcdt_min_angle;
  double rcdt_max_edge;
  Matrix<double> quality;
  int rcdt_max_n0;
  int rcdt_max_n1;
  bool rcdt;

public:
  Options(Rcpp::List& options, size_t rows) :
  cutoff(1.0e-12),
  sphere_tolerance(1.0e-7),
  cet_sides(8),
  cet_margin(-0.1),
  rcdt_min_angle(21),
  rcdt_max_edge(-1.0),
  quality(1),
  rcdt_max_n0(-1),
  rcdt_max_n1(-1),
  rcdt(true) {
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "cutoff"))
      cutoff = options["cutoff"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "sphere_tolerance"))
      sphere_tolerance = options["sphere_tolerance"];
    if (Rcpp_is_element<Rcpp::IntegerVector>(options, "cet_sides"))
      cet_sides = options["cet_sides"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "cet_margin"))
      cet_margin = options["cet_margin"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "rcdt_min_angle"))
      rcdt_min_angle = options["rcdt_min_angle"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "rcdt_max_edge"))
      rcdt_max_edge = options["rcdt_max_edge"];
    if (Rcpp_is_element<Rcpp::LogicalVector>(options, "rcdt"))
      rcdt = options["rcdt"];

    /* Construct quality info */
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "quality")) {
      quality = Rcpp::as<Rcpp::NumericVector>(options["quality"]);
      for (size_t r = quality.rows(); r < rows; r++)
        quality(r, 0) = rcdt_max_edge;
      quality.rows(rows); /* Make sure we have the right number of rows */
    } else {
      quality.rows(rows);
      for (size_t r = 0; r < rows; r++)
        quality(r, 0) = rcdt_max_edge;
    }

    if (Rcpp_is_element<Rcpp::IntegerVector>(options, "rcdt_max_n0"))
      rcdt_max_n0 = options["rcdt_max_n0"];
    if (Rcpp_is_element<Rcpp::IntegerVector>(options, "rcdt_max_n1"))
      rcdt_max_n1 = options["rcdt_max_n1"];
  };

};



Mesh Rcpp_import_mesh(Rcpp::NumericMatrix mesh_loc,
                      Rcpp::IntegerMatrix mesh_tv,
                      MatrixC & matrices,
                      Rcpp::List options) {
  const bool useVT = true;
  const bool useTTi = true;

  matrices.attach("mesh_loc", new Matrix3double(Matrix<double>(mesh_loc)), true);
  FMLOG("'mesh_loc' points imported." << std::endl);
  matrices.attach("mesh_tv", new Matrix<int>(mesh_tv), true);
  FMLOG("'mesh_tv' points imported." << std::endl);

  Matrix<double>& iS0 = matrices.DD("mesh_loc");
  Matrix<int>& TV0 = matrices.DI("mesh_tv");

  /* Initialise mesh structure */
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  //  if ((iS0.rows() > 0) && (iS0.cols() < 2)) {
  //    /* 1D data. Not implemented */
  //    FMLOG("1D data not implemented." << std::endl);
  //    return Rcpp::List();
  //  }

  if (iS0.rows() > 0) {
    //    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(iS0);
  }

  Options the_options(options, iS0.rows());

  //  double sphere_tolerance = 1e-10;
  (void)M.auto_type(the_options.sphere_tolerance);

  M.TV_set(TV0);

  return M;
}







// #include "qtool.h"
//
// //' @title Compute sparse matrix inverse
// //'
// //' @description
// //' Requires RcppEigen which is not compiled in by default
// //'
// //' @param AA A sparse matrix
// //' @keywords internal
// // [[Rcpp::export]]
// Rcpp::List C_qinv(SEXP AA) {
// #ifdef FMESHER_WITH_EIGEN
//   const EigenMSM<double> A(Rcpp::as<EigenMSM<double>>(AA));
//
//   QTool<double> Q;
//   Q.Q(A);
//
//   Rcpp::List ret;
//   ret["Qinv"] = Q.S();
//   return ret;
// #else
//   Rcpp::stop("Unsupported method C_qinv; fmesher was built without FMESHER_WITH_EIGEN");
//   Rcpp::List ret;
//   return ret;
// #endif
// }



//' @title Globe points
//'
//' @description
//' Create points on a globe
//'
//' @param globe integer; the number of edge subdivision segments, 1 or higher
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix fmesher_globe_points(Rcpp::IntegerVector globe) {
  MatrixC matrices;

  matrices.attach(
    ".globe",
    (Matrix<double> *)fmesh::make_globe_points(globe[0], 1.0),
    true);
  FMLOG("globe points constructed." << std::endl);

  return Rcpp::wrap(matrices.DD(".globe"));
}




//' @title Refined Constrained Delaunay Triangulation
//'
//' @description
//' (...)
//'
//' @param options list of triangulation options
//' @param loc numeric matrix; initial points to include
//' @param tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param boundary 2-column integer matrix with 0-based vertex indices for each
//' boundary edge constraint
//' @param interior 2-column integer matrix with 0-based vertex indices for each
//' interior edge constraint
//' @param boundary_grp integer vector with group labels
//' @param interior_grp integer vector with group labels
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_rcdt(Rcpp::List options,
                        Rcpp::NumericMatrix loc,
                        Rcpp::Nullable<Rcpp::IntegerMatrix> tv = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerMatrix> boundary = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerMatrix> interior = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerVector> boundary_grp = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerVector> interior_grp = R_NilValue) {
  const bool useVT = true;
  const bool useTTi = true;

  MatrixC matrices;

  matrices.attach("loc", new Matrix<double>(loc), true);
  FMLOG("'loc' points imported." << std::endl);

  Matrix<double>& iS0 = matrices.DD("loc");
  Matrix<int>* TV0 = NULL;
  if (!tv.isNull()) {
    matrices.attach("tv0",
                    new Matrix<int>(Rcpp::as<Rcpp::IntegerMatrix>(tv)),
                    true);
    FMLOG("'tv0' points imported." << std::endl);
    TV0 = &matrices.DI("tv0");
  }

  Options rcdt_options(options, iS0.rows());
  FMLOG("rcdt_options parsed" << std::endl);

  /* Prepare boundary/interior edges */
  matrices.attach("boundary", new Matrix<int>(2), true);
  matrices.attach("interior", new Matrix<int>(2), true);
  matrices.attach("boundary_grp", new Matrix<int>(1), true);
  matrices.attach("interior_grp", new Matrix<int>(1), true);
  if (!boundary.isNull()) {
    matrices.DI("boundary") = Rcpp::as<Rcpp::IntegerMatrix>(boundary);
  }
  if (!interior.isNull()) {
    matrices.DI("interior") = Rcpp::as<Rcpp::IntegerMatrix>(interior);
  }
  if (!boundary_grp.isNull()) {
    matrices.DI("boundary_grp") = Rcpp::as<Rcpp::IntegerVector>(boundary_grp);
  } else {
    matrices.DI("boundary_grp")(0, 0) = 1;
  }
  if (!interior_grp.isNull()) {
    matrices.DI("interior_grp") = Rcpp::as<Rcpp::IntegerVector>(interior_grp);
  } else {
    matrices.DI("interior_grp")(0, 0) = 1;
  }

  constrListT cdt_boundary;
  constrListT cdt_interior;
  if (!boundary.isNull()) {
    prepare_cdt_input(matrices.DI("boundary"),
                      matrices.DI("boundary_grp"),
                      cdt_boundary);
  }
  if (!interior.isNull()) {
    prepare_cdt_input(matrices.DI("interior"),
                      matrices.DI("interior_grp"),
                      cdt_interior);
  }

  /* Prepare to filter out points at distance not greater than 'cutoff' */
  matrices.attach("idx", new Matrix<int>(iS0.rows(), 1), true);
  matrices.output("idx");
  Matrix<int> &idx = matrices.DI("idx").clear();

  filter_locations(iS0, idx, rcdt_options.cutoff);

  /* Remap vertex input references */
  if (TV0) {
    remap_vertex_indices(idx, *TV0);
  }
  remap_vertex_indices(idx, cdt_boundary);
  remap_vertex_indices(idx, cdt_interior);

  /* Initialise mesh structure */
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  if ((iS0.rows() > 0) && (iS0.cols() < 2)) {
    /* 1D data. Not implemented */
    FMLOG("1D data not implemented." << std::endl);
    return Rcpp::List();
  }

  if (iS0.rows() > 0) {
    FMLOG("Append S0." << std::endl);
    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(S0);
  }

  FMLOG("Auto-detect manifold type." << std::endl);
  //  double sphere_tolerance = 1e-10;
  (void)M.auto_type(rcdt_options.sphere_tolerance);

  if (TV0) {
    FMLOG("Set TV0." << std::endl);
    M.TV_set(*TV0);
  }

  FMLOG("Attach 's'." << std::endl);
  matrices.attach(string("s"), &M.S(), false);
  FMLOG("Attach 'tv'." << std::endl);
  matrices.attach("tv", &M.TV(), false);
  FMLOG("Set output of 's' and 'tv'." << std::endl);
  matrices.output("s").output("tv");

  FMLOG("Create MeshC helper." << std::endl);
  MeshC MC(&M);
  MC.setOptions(MC.getOptions() | MeshC::Option_offcenter_steiner);

  if ((M.type() != Mesh::Mtype_plane) &&
      (M.type() != Mesh::Mtype_sphere)) {
    if (M.nT() == 0) {
      FMLOG_(
        "Points not in the plane or on a sphere, and triangulation empty."
        << std::endl);
    }
    /* Remove everything outside the boundary segments, if any. */
    MC.PruneExterior();
    invalidate_unused_vertex_indices(M, idx);
    /* Nothing more to do here.  Cannot refine non R2/S2 meshes. */
  } else {
    /* If we don't already have a triangulation, we must create one. */
    if (M.nT() == 0) {
      FMLOG("cet_sides = " << rcdt_options.cet_sides << std::endl);
      FMLOG("cet_margin = " << rcdt_options.cet_margin << std::endl);
      if (!MC.CET(rcdt_options.cet_sides, rcdt_options.cet_margin)) {
        FMLOG_("CET creation failed, exiting." << std::endl);
        return Rcpp::wrap(matrices);
      }
    }

    /* It is more robust to add the constraints before the rest of the
     nodes are added.  This allows points to fall onto constraint
     segments, subdividing them as needed. */
    if (cdt_boundary.size() > 0)
      MC.CDTBoundary(cdt_boundary);
    if (cdt_interior.size() > 0)
      MC.CDTInterior(cdt_interior);

    /* Add the rest of the nodes. */
    vertexListT vertices;
    for (size_t v = 0; v < iS0.rows(); v++)
      vertices.push_back(v);

    MC.DT(vertices);

    /* Remove everything outside the boundary segments, if any. */
    MC.PruneExterior();
    invalidate_unused_vertex_indices(M, idx);

    if ((rcdt_options.rcdt) &&
        (rcdt_options.rcdt_max_edge > 0)) {
      /* Calculate the RCDT: */
      MC.RCDT(rcdt_options.rcdt_min_angle,
              rcdt_options.rcdt_max_edge,
              rcdt_options.quality.raw(),
              rcdt_options.quality.rows(),
              rcdt_options.rcdt_max_n0,
              rcdt_options.rcdt_max_n1);
      FMLOG(MC << endl);
    }
    /* Done constructing the triangulation. */

    /* Calculate and collect output. */

    matrices.attach("segm.bnd.idx", new Matrix<int>(2), true,
                    fmesh::IOMatrixtype_general);
    matrices.attach("segm.bnd.grp", new Matrix<int>(1), true,
                    fmesh::IOMatrixtype_general);
    MC.segments(true,
                &matrices.DI("segm.bnd.idx"),
                &matrices.DI("segm.bnd.grp"));

    matrices.output("segm.bnd.idx").output("segm.bnd.grp");

    matrices.attach("segm.int.idx", new Matrix<int>(2), true,
                    fmesh::IOMatrixtype_general);
    matrices.attach("segm.int.grp", new Matrix<int>(1), true,
                    fmesh::IOMatrixtype_general);
    MC.segments(false, &matrices.DI("segm.int.idx"),
                &matrices.DI("segm.int.grp"));

    matrices.output("segm.int.idx").output("segm.int.grp");
  }

  matrices.attach("tt", &M.TT(), false);
  M.useVT(true);
  matrices.attach("vt", &M.VT(), false);
  M.useTTi(true);
  matrices.attach("tti", &M.TTi(), false);
  matrices.attach("vv", new SparseMatrix<int>(M.VV()), true,
                  fmesh::IOMatrixtype_symmetric);

  matrices.output("tt").output("tti").output("vt").output("vv");

//  FMLOG("Manifold output." << std::endl);
//  /* Output the manifold type. */
//  matrices.attach("manifold", new Matrix<int>(1), true,
//                  fmesh::IOMatrixtype_general);
//  Matrix<int> &manifold = matrices.DI("manifold");
//  manifold(0, 0) = M.type();
//  matrices.output("manifold");

  Rcpp::List out = Rcpp::wrap(matrices);

  switch (M.type()) {
  case Mesh::Mtype_manifold:
    out["manifold"] = "M2";
    break;
  case Mesh::Mtype_plane:
    out["manifold"] = "R2";
    break;
  case Mesh::Mtype_sphere:
    out["manifold"] = "S2";
    break;
  }

  return out;
}




//' @title Barycentric coordinate computation
//'
//' @description
//' Locate points and compute triangular barycentric coordinates
//'
//' @param loc numeric matrix; coordinates of points to locate in the mesh
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param options list of triangulation options
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' b <- fmesher_bary(m$s,
//'                   m$tv,
//'                   matrix(c(0.5, 0.5), 1, 2),
//'                   list())
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_bary(Rcpp::NumericMatrix mesh_loc,
                        Rcpp::IntegerMatrix mesh_tv,
                        Rcpp::NumericMatrix loc,
                        Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, options);
  Options rcdt_options(options, M.nV());

  FMLOG("barycentric coordinate output." << std::endl);
  if ((M.type() != Mesh::Mtype_plane) &&
      (M.type() != Mesh::Mtype_sphere)) {
    FMLOG_("Cannot calculate points2mesh mapping for non R2/S2 manifolds"
             << std::endl);
    return Rcpp::List();
  }

  matrices.attach("loc", new Matrix3double(Matrix<double>(loc)), true);
  Matrix<double>& points2mesh = matrices.DD("loc");

  size_t points_n = points2mesh.rows();
  Matrix<int> &points2mesh_t =
    matrices.attach(string("t"), new Matrix<int>(points_n, 1), true);
  Matrix<double> &points2mesh_b = matrices.attach(
    string("bary"), new Matrix<double>(points_n, 3), true);
  matrices.matrixtype("t", fmesh::IOMatrixtype_general);
  matrices.matrixtype("bary", fmesh::IOMatrixtype_general);
  matrices.output("t").output("bary");

  map_points_to_mesh(M, points2mesh, points2mesh_t, points2mesh_b);

  return Rcpp::wrap(matrices);
}




//' @title Rotationally invariant spherical B-splines
//'
//' @description
//' Compute rotationally invariant spherical B-splines on the unit sphere
//'
//' @param loc numeric vector/matrix; coordinates of points to locate in the mesh,
//' only the z-coordinates are used (`sin(latitude)`)
//' @param n The number of basis functions
//' @param degree The polynomial basis degree
//' @param uniform logical; If `TRUE`, the knots are spaced uniformly by latitude,
//' if `FALSE`, the knots are spaced uniformly by `sin(latitude)`
//' @rdname fmesher_spherical_bsplines
//' @examples
//' m <- fm_rcdt_2d(globe = 1)
//' fmesher_spherical_bsplines(m$loc, n = 3, degree = 2, uniform = FALSE)
//' fmesher_spherical_bsplines1(m$loc[, 3], n = 3, degree = 2, uniform = FALSE)
//' @export
//' @keywords internal
// [[Rcpp::export]]
SEXP fmesher_spherical_bsplines1(Rcpp::NumericVector loc,
                                 int n,
                                 int degree,
                                 Rcpp::LogicalVector uniform) {
  if (n < 0) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (degree < 1) {
    Rcpp::stop("'degree' must be at least 0.");
  }
  if (n <= degree) {
    Rcpp::stop("'n' must be larger than 'degree'");
  }

  MatrixC matrices;
  matrices.attach("loc", new Matrix<double>(loc), true);

  FMLOG("bspline output." << std::endl);

  bool bool_uniform = Rcpp::is_true(Rcpp::all(uniform));
  if (bool_uniform) {
    FMLOG("uniform = TRUE" << std::endl);
  } else {
    FMLOG("uniform = FALSE" << std::endl);
  }
  matrices.attach(
    string("bspline"),
    spherical_bsplines1(matrices.DD("loc"), n, degree, bool_uniform),
    true);
  matrices.matrixtype("bspline", fmesh::IOMatrixtype_general);
  matrices.output("bspline");

  return Rcpp::wrap(matrices.DD("bspline"));
}

//' @rdname fmesher_spherical_bsplines
//' @export
// [[Rcpp::export]]
SEXP fmesher_spherical_bsplines(Rcpp::NumericMatrix loc,
                                int n,
                                int degree,
                                Rcpp::LogicalVector uniform) {
  if (n < 0) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (degree < 1) {
    Rcpp::stop("'degree' must be at least 0.");
  }
  if (n <= degree) {
    Rcpp::stop("'n' must be larger than 'degree'");
  }
  if (loc.cols() < 3) {
    Rcpp::stop("'ncol(loc)' must be at least 3.");
  }

  MatrixC matrices;
  matrices.attach("loc", new Matrix3double(Matrix<double>(loc)), true);

  FMLOG("bspline output." << std::endl);

  bool bool_uniform = Rcpp::is_true(Rcpp::all(uniform));
  if (bool_uniform) {
    FMLOG("uniform = TRUE" << std::endl);
  } else {
    FMLOG("uniform = FALSE" << std::endl);
  }
  matrices.attach(
    string("bspline"),
    spherical_bsplines(matrices.DD("loc"), n, degree, bool_uniform),
    true);
  matrices.matrixtype("bspline", fmesh::IOMatrixtype_general);
  matrices.output("bspline");

  return Rcpp::wrap(matrices.DD("bspline"));
  //  return Rcpp::wrap(matrices);
}




//' @title Finite element matrix computation
//'
//' @description
//' Construct finite element structure matrices
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param fem_order_max integer; the highest operator order to compute
//' @param aniso If non-NULL, a `list(gamma, v)`. Calculates anisotropic structure
//' matrices (in addition to the regular) for \eqn{\gamma}{gamma} and \eqn{v}{v} for
//' an anisotropic operator \eqn{\nabla\cdot H \nabla}{div H grad}, where
//' \eqn{H=\gamma I + v v^\top}{H = gamma I + v v'}.
//' Currently (2023-08-05) the fields need to be given per vertex.
//' @param options list of triangulation options (`sphere_tolerance`)
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' b <- fmesher_fem(m$s, m$tv, fem_order_max = 2, aniso = NULL, options = list())
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_fem(Rcpp::NumericMatrix mesh_loc,
                       Rcpp::IntegerMatrix mesh_tv,
                       int fem_order_max,
                       Rcpp::Nullable<Rcpp::List> aniso,
                       Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, options);

  FMLOG("Compute finite element matrices." << std::endl);

  if (fem_order_max >= 0) {
    FMLOG("fem output." << std::endl)
    SparseMatrix<double> &C0 = matrices.SD("c0").clear();
    SparseMatrix<double> &C1 = matrices.SD("c1").clear();
    SparseMatrix<double> &B1 = matrices.SD("b1").clear();
    SparseMatrix<double> &G = matrices.SD("g1").clear();
    SparseMatrix<double> &K = matrices.SD("k1").clear();
    /* K1=G1-B1, K2=K1*inv(C0)*K1, ... */
    Matrix<double> &Tareas = matrices.DD("ta").clear();

    M.calcQblocks(C0, C1, G, B1, Tareas);

    matrices.attach(string("va"), new Matrix<double>(diag(C0)), true);

    K = G - B1;

    matrices.matrixtype("c0", fmesh::IOMatrixtype_diagonal);
    matrices.matrixtype("c1", fmesh::IOMatrixtype_symmetric);
    matrices.matrixtype("b1", fmesh::IOMatrixtype_general);
    matrices.matrixtype("g1", fmesh::IOMatrixtype_symmetric);
    matrices.matrixtype("k1", fmesh::IOMatrixtype_general);
    matrices.output("c0");
    matrices.output("c1");
    matrices.output("b1");
    matrices.output("g1");
    matrices.output("k1");
    matrices.output("va");
    matrices.output("ta");

    SparseMatrix<double> C0inv = inverse(C0, true);
    // Protect temporary local variables
    {
      SparseMatrix<double> tmp = G * C0inv;
      SparseMatrix<double> *a;
      SparseMatrix<double> *b = &G;
      for (size_t i = 1; int(i) < fem_order_max; i++) {
        std::stringstream ss;
        ss << i + 1;
        std::string Gname = "g" + ss.str();
        a = b;
        b = &(matrices.SD(Gname).clear());
        *b = tmp * (*a);
        matrices.matrixtype(Gname, fmesh::IOMatrixtype_symmetric);
        matrices.output(Gname);
      }
      tmp = C0inv * K;
      b = &K;
      for (size_t i = 1; int(i) < fem_order_max; i++) {
        std::stringstream ss;
        ss << i + 1;
        std::string Kname = "k" + ss.str();
        a = b;
        b = &(matrices.SD(Kname).clear());
        *b = (*a) * tmp;
        matrices.matrixtype(Kname, fmesh::IOMatrixtype_general);
        matrices.output(Kname);
      }
    }

    if (!aniso.isNull()) {
      FMLOG("Compute anisotropic finite element matrices." << std::endl);
      if (Rcpp::as<Rcpp::List>(aniso).size() < 2) {
        Rcpp::stop("'aniso' list must have at least two elements.");
      }
      matrices.attach("gamma_field",
                      new Matrix<double>(
                          Rcpp::as<Rcpp::NumericVector>(
                            Rcpp::as<Rcpp::List>(aniso)[0]
                          )
                      ),
                      true);
      FMLOG("'gamma_field' imported." << std::endl);
      matrices.attach("vector_field",
                      new Matrix<double>(
                          Rcpp::as<Rcpp::NumericMatrix>(
                            Rcpp::as<Rcpp::List>(aniso)[1]
                          )
                      ),
                      true);
      FMLOG("'vector_field' imported." << std::endl);
      if (matrices.DD("gamma_field").rows() < M.nV()) {
        Rcpp::stop("'aniso[[1]]' length should match the number of vertices.");
      }
      if (matrices.DD("vector_field").rows() < M.nV()) {
        Rcpp::stop("'aniso[[2]]' rows should match the number of vertices.");
      }

      SparseMatrix<double> &Gani = matrices.SD("g1aniso").clear();
      M.calcQblocksAni(Gani,
                       matrices.DD("gamma_field"),
                       matrices.DD("vector_field"));
      matrices.output("g1aniso");

      // Protect temporary local variables
      {
        SparseMatrix<double> tmp = Gani * C0inv;
        SparseMatrix<double> *a;
        SparseMatrix<double> *b = &Gani;
        for (size_t i = 1; int(i) < fem_order_max; i++) {
          std::stringstream ss;
          ss << i + 1;
          std::string Gname = "g" + ss.str() + "aniso";
          a = b;
          b = &(matrices.SD(Gname).clear());
          *b = tmp * (*a);
          matrices.matrixtype(Gname, fmesh::IOMatrixtype_symmetric);
          matrices.output(Gname);
        }
      }
    }

  }

  return Rcpp::wrap(matrices);
}




//' @title Split lines at triangle edges
//'
//' @description
//' (...)
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param loc numeric coordinate matrix
//' @param idx 2-column integer matrix
//' @param options list of triangulation options (`sphere_tolerance`)
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_split_lines(
    Rcpp::NumericMatrix mesh_loc,
    Rcpp::IntegerMatrix mesh_tv,
    Rcpp::NumericMatrix loc,
    Rcpp::IntegerMatrix idx,
    Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, options);

  FMLOG("Compute line splitting." << std::endl);

  matrices.attach("loc", new Matrix3double(Matrix<double>(loc)), true);
  matrices.attach("idx", new Matrix<int>(idx), true);

  /* Make sure we have a Nx3 matrix: */
  Matrix<double> *splitloc1 = new Matrix<double>(3);
  Matrix<int> *splitidx1 = new Matrix<int>(2);
  Matrix<int> *splittriangle1 = new Matrix<int>(1);
  Matrix<double> *splitbary1 = new Matrix<double>(3);
  Matrix<double> *splitbary2 = new Matrix<double>(3);
  Matrix<int> *splitorigin1 = new Matrix<int>(1);

  split_line_segments_on_triangles(
    M, matrices.DD("loc"), matrices.DI("idx"), *splitloc1,
    *splitidx1, *splittriangle1, *splitbary1, *splitbary2, *splitorigin1);

  /* Now it's ok to overwrite potential input split* matrices. */
  matrices.attach("split.loc", splitloc1, true);
  matrices.attach("split.idx", splitidx1, true);
  matrices.attach("split.t", splittriangle1, true);
  matrices.attach("split.b1", splitbary1, true);
  matrices.attach("split.b2", splitbary2, true);
  matrices.attach("split.origin", splitorigin1, true);
  matrices.output("split.loc").output("split.idx");
  matrices.output("split.b1").output("split.b2");
  matrices.output("split.t").output("split.origin");

  return Rcpp::wrap(matrices);
}



// //' @title Test the matrix I/O system
// //'
// //' @param args_input Input argument list
// //' @examples
// //' \dontrun{
// //' A <- Matrix::sparseMatrix(i=1:4,j=4:1,x=2:5,dims=c(4,4))
// //' inp <- list(
// //'   A = fm_as_dgTMatrix(A),
// //'   Bd = matrix((11:22)+0.5,4,3),
// //'   Bi = matrix(121L:132L,4,3),
// //'   B1d=as.matrix((31:34)+0.5),
// //'   B1i=as.matrix(41L:44L),
// //'   Ad = fm_as_fmesher_sparse(A)
// //' )
// //' inp[["BdM"]] <- as(inp[["Bd"]], "unpackedMatrix")
// //' out <- C_matrixio_test2(args_input = inp)
// //' str(out)
// //' }
// //' @keywords internal
// // [[Rcpp::export]]
// Rcpp::List C_matrixio_test2(Rcpp::List args_input) {
//   MatrixC matrices(args_input);
//   Rcpp::List ret = Rcpp::wrap(matrices.output("-"));
//   return (ret);
// }



// //' @title Test the matrix I/O system
// //'
// //' @param args_input Input argument list
// //' @examples
// //' \dontrun{
// //' A <- Matrix::sparseMatrix(i=1:4,j=4:1,x=2:5,dims=c(4,4))
// //' out <- C_matrixio_test(args_input=list(
// //'   A = fm_as_dgTMatrix(A),
// //'   Bd = matrix((11:22)+0.5,4,3),
// //'   Bi = matrix(121L:132L,4,3),
// //'   B1d=as.matrix((31:34)+0.5),
// //'   B1i=as.matrix(41L:44L),
// //'   Ad = fm_as_fmesher_sparse(A)
// //' ))
// //' Aout <- fm_as_dgTMatrix(out[["Ad"]])
// //' A
// //' Aout
// //' }
// //' @keywords internal
// // [[Rcpp::export]]
//      Rcpp::List C_matrixio_test(Rcpp::List args_input) {
//        MatrixC matrices;
//
//        //  matrices.attach("loc", new Matrix<double>(Rcpp::as<EigenMM<double>>(args_input["loc"])), true);
//        //  matrices.attach("tv", new Matrix<int>(Rcpp::as<EigenMM<int>>(args_input["tv"])), true);
//
//        bool is_list = Rcpp::is<Rcpp::List>(args_input);
//        bool is_numeric_matrix = Rcpp::is<Rcpp::NumericMatrix>(args_input["A"]);
//        bool is_numeric_vector = Rcpp::is<Rcpp::NumericVector>(args_input["A"]);
//        bool is_integer_matrix = Rcpp::is<Rcpp::IntegerMatrix>(args_input["A"]);
//        bool is_integer_vector = Rcpp::is<Rcpp::IntegerVector>(args_input["A"]);
//
//        Rcpp::NumericMatrix Bd = Rcpp::as<Rcpp::NumericMatrix>(args_input["Bd"]);
//        Rcpp::IntegerMatrix Bi = Rcpp::as<Rcpp::IntegerMatrix>(args_input["Bi"]);
//        Rcpp::NumericVector B1d = Rcpp::as<Rcpp::NumericVector>(args_input["B1d"]);
//        Rcpp::IntegerVector B1i = Rcpp::as<Rcpp::IntegerVector>(args_input["B1i"]);
//
//
//        fmesh::Matrix<double> Bdd = Bd;
//        fmesh::Matrix<int> Bdi(Rcpp::as<Rcpp::IntegerMatrix>(Bd));
//        fmesh::Matrix<double> Bid(Rcpp::as<Rcpp::NumericMatrix>(Bi));
//        fmesh::Matrix<int> Bii(Bi);
//
//        FMLOG_("Bdd: " << Bdd << std::endl);
//        FMLOG_("Bdi: " << Bdi << std::endl);
//        FMLOG_("Bid: " << Bid << std::endl);
//        FMLOG_("Bii: " << Bii << std::endl);
//
//        fmesh::Matrix1<double> Bdd1 = B1d;
//        fmesh::Matrix1<double> Bdd_1 = Rcpp::NumericVector(Bd(Rcpp::_, 1));
//        fmesh::Matrix1<int> Bdi1(Rcpp::as<Rcpp::IntegerVector>(B1d));
//        fmesh::Matrix1<double> Bid1(Rcpp::as<Rcpp::NumericVector>(B1i));
//        fmesh::Matrix1<int> Bii1(B1i);
//
//        FMLOG_("Bdd1: " << Bdd1 << std::endl);
//        FMLOG_("Bdd_1: " << Bdd_1 << std::endl);
//        FMLOG_("Bdi1: " << Bdi1 << std::endl);
//        FMLOG_("Bid1: " << Bid1 << std::endl);
//        FMLOG_("Bii1: " << Bii1 << std::endl);
//
//        fmesh::Matrix3<double> Bdd3 = Bd;
//        fmesh::Matrix3<int> Bdi3(Rcpp::as<Rcpp::IntegerMatrix>(Bd));
//        fmesh::Matrix3<double> Bid3(Rcpp::as<Rcpp::NumericMatrix>(Bi));
//        fmesh::Matrix3<int> Bii3(Bi);
//
//        FMLOG_("Bdd3: " << Bdd3 << std::endl);
//        FMLOG_("Bdi3: " << Bdi3 << std::endl);
//        FMLOG_("Bid3: " << Bid3 << std::endl);
//        FMLOG_("Bii3: " << Bii3 << std::endl);
//
//        const Rcpp::List Ad(Rcpp::as<Rcpp::List>(args_input["Ad"]));
//
//        fmesh::SparseMatrix<double> Ad_fm(Ad);
//
//        //  const EigenMSM<int> Ai(Rcpp::as<EigenMSM<int>>(args_input["Ai"]));
//
//        //  bool is_msm = Rcpp::is<Eigen::SparseMatrix<double>>(args_input["a"]);
//
//        matrices.attach("Ad_fm", &Ad_fm, false);
//        matrices.output("Ad_fm");
//
//        MatrixC mat2(args_input);
//
//        FMLOG_(mat2.DI("tv"))
//
//          Rcpp::List ret;
//        ret["is_list"] = is_list;
//        ret["is_numeric_matrix"] = is_numeric_matrix;
//        ret["is_numeric_vector"] = is_numeric_vector;
//        ret["is_integer_matrix"] = is_integer_matrix;
//        ret["is_integer_vector"] = is_integer_vector;
//        //  ret["A"] = A;
//        ret["Ad"] = Ad;
// #ifdef FMESHER_WITH_EIGEN
//        ret["Ad_fm"] = Ad_fm.EigenSparseMatrix();
// #endif
//        ret["Ad_fm_auto"] = Ad_fm;
//        ret["Ad_fm_ijx"] = Ad_fm.fmesher_sparse();
//        ret["Bid3"] = Bid3;
//        ret["Bdi3"] = Bdi3;
//        ret["matrices"] = matrices;
//        ret["mat2"] = mat2.output("-");
//        return (ret);
//      }




#endif
