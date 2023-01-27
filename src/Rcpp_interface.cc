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
using fmesh::MatrixC;
using fmesh::Mesh;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::SparseMatrix;
using fmesh::TriangleLocator;
using fmesh::Vector3;
using fmesh::vertexListT;

template <class T> using EigenM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <class T> using EigenM1 = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <class T> using EigenSM = Eigen::SparseMatrix<T>;
template <class T> using EigenMM = Eigen::Map<EigenM<T>>;
template <class T> using EigenMM1 = Eigen::Map<EigenM1<T>>;
template <class T> using EigenMSM = Eigen::Map<EigenSM<T>>;

// const bool useVT = true;
// const bool useTTi = true;

#include "qtool.h"

//' Compute sparse matrix inverse
//'
//' @param AA A sparse matrix
// [[Rcpp::export]]
Rcpp::List C_qinv(SEXP AA) {
  // Eigen::SparseMatrix<double> C_qinv(SEXP AA)
  const EigenMSM<double> A(Rcpp::as<EigenMSM<double>>(AA));

  QTool<double> Q;
  Q.Q(A);

  Rcpp::List ret;
  ret["Qinv"] = Q.S();
  return (ret);
  //  return Rcpp::List::create(Rcpp::Named("Q") = Q.S());
  //  return Q.S();
}

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




template <typename T>
bool is_element(const Rcpp::List& list, std::string name) {
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
    if (is_element<Rcpp::NumericVector>(options, "cutoff"))
      cutoff = options["cutoff"];
    if (is_element<Rcpp::NumericVector>(options, "sphere_tolerance"))
      sphere_tolerance = options["sphere_tolerance"];
    if (is_element<Rcpp::IntegerVector>(options, "cet_sides"))
      cet_sides = options["cet_sides"];
    if (is_element<Rcpp::NumericVector>(options, "cet_margin"))
      cet_margin = options["cet_margin"];
    if (is_element<Rcpp::NumericVector>(options, "rcdt_min_angle"))
      rcdt_min_angle = options["rcdt_min_angle"];
    if (is_element<Rcpp::NumericVector>(options, "rcdt_max_edge"))
      rcdt_max_edge = options["rcdt_max_edge"];
    if (is_element<Rcpp::LogicalVector>(options, "rcdt"))
      rcdt = options["rcdt"];

    /* Construct quality info */
    if (is_element<Rcpp::NumericVector>(options, "quality")) {
      quality = Rcpp::as<Rcpp::NumericVector>(options["quality"]);
      for (int r = quality.rows(); r < rows; r++)
        quality(r, 0) = rcdt_max_edge;
      quality.rows(rows); /* Make sure we have the right number of rows */
    } else {
      quality.rows(rows);
      for (int r = 0; r < rows; r++)
        quality(r, 0) = rcdt_max_edge;
    }

    if (is_element<Rcpp::IntegerVector>(options, "rcdt_max_n0"))
      rcdt_max_n0 = options["rcdt_max_n0"];
    if (is_element<Rcpp::IntegerVector>(options, "rcdt_max_n1"))
      rcdt_max_n1 = options["rcdt_max_n1"];
  };

};



//' @title Refined Constrained Delaunay Triangulation
//'
//' @description
//' (...)
//'
//' @param options list of triangulation options
//' @param loc numeric matrix; initial points to include
//' @param tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param boundary 2-column integer matrix with 0-based vertex indices for each
//' boundary edge constrain
//' @param interior 2-column integer matrix with 0-based vertex indices for each
//' interior edge constraint
//' @param boundary_grp integer vector with group lables
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
                    new Matrix<double>(Rcpp::as<Rcpp::IntegerMatrix>(tv)),
                    true);
    FMLOG("'tv0' points imported." << std::endl);
    TV0 = &matrices.DI("tv0");
  }

  Options rcdt_options(options, iS0.rows());

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
    matrices.DI("boundary_grp") = Rcpp::as<Rcpp::IntegerMatrix>(boundary_grp);
  } else {
    matrices.DI("boundary_grp")(0, 0) = 1;
  }
  if (!interior_grp.isNull()) {
    matrices.DI("interior_grp") = Rcpp::as<Rcpp::IntegerMatrix>(interior_grp);
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
    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(S0);
  }

  //  double sphere_tolerance = 1e-10;
  (void)M.auto_type(rcdt_options.sphere_tolerance);

  if (TV0) {
    M.TV_set(*TV0);
  }

  matrices.attach(string("s"), &M.S(), false);
  matrices.attach("tv", &M.TV(), false);
  matrices.output("s").output("tv");

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
      MC.CET(rcdt_options.cet_sides, rcdt_options.cet_margin);
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
//' (...)
//'
//' @param loc numeric matrix; coordinates of points to locate in the mesh
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param options list of triangulation options
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' b <- fmesher_bary(matrix(c(0.5, 0.5), 1, 2),
//'                   m$s,
//'                   m$tv,
//'                   list())
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_bary(Rcpp::NumericMatrix loc,
                        Rcpp::NumericMatrix mesh_loc,
                        Rcpp::IntegerMatrix mesh_tv,
                        Rcpp::List options) {
  const bool useVT = true;
  const bool useTTi = true;

  MatrixC matrices;

  matrices.attach("loc", new Matrix<double>(loc), true);
  FMLOG("'loc' points imported." << std::endl);
  matrices.attach("mesh_loc", new Matrix<double>(mesh_loc), true);
  FMLOG("'mesh_loc' points imported." << std::endl);
  matrices.attach("mesh_tv", new Matrix<int>(mesh_tv), true);
  FMLOG("'mesh_tv' points imported." << std::endl);

  Matrix<double>& points2mesh = matrices.DD("loc");
  Matrix<double>& iS0 = matrices.DD("mesh_loc");
  Matrix<int>& TV0 = matrices.DI("mesh_tv");

  Options rcdt_options(options, iS0.rows());

  /* Initialise mesh structure */
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  if ((iS0.rows() > 0) && (iS0.cols() < 2)) {
    /* 1D data. Not implemented */
    FMLOG("1D data not implemented." << std::endl);
    return Rcpp::List();
  }

  if (iS0.rows() > 0) {
    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(S0);
  }

  //  double sphere_tolerance = 1e-10;
  (void)M.auto_type(rcdt_options.sphere_tolerance);

  M.TV_set(TV0);

  FMLOG("barycentric coordinate output." << std::endl);
  if ((M.type() != Mesh::Mtype_plane) &&
      (M.type() != Mesh::Mtype_sphere)) {
    FMLOG_("Cannot calculate points2mesh mapping for non R2/S2 manifolds"
             << std::endl);
    return Rcpp::List();
  }

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




//' Test the matrix I/O system
//'
//' @param args_input Input argument list
// [[Rcpp::export]]
Rcpp::List C_matrixio_test(Rcpp::List args_input) {
  // Eigen::SparseMatrix<double> C_qinv(SEXP AA)

  MatrixC matrices;

//  matrices.attach("loc", new Matrix<double>(Rcpp::as<EigenMM<double>>(args_input["loc"])), true);
//  matrices.attach("tv", new Matrix<int>(Rcpp::as<EigenMM<int>>(args_input["tv"])), true);

  bool is_list = Rcpp::is<Rcpp::List>(args_input);
  bool is_numeric_matrix = Rcpp::is<Rcpp::NumericMatrix>(args_input["A"]);
  bool is_numeric_vector = Rcpp::is<Rcpp::NumericVector>(args_input["A"]);
  bool is_integer_matrix = Rcpp::is<Rcpp::IntegerMatrix>(args_input["A"]);
  bool is_integer_vector = Rcpp::is<Rcpp::IntegerVector>(args_input["A"]);

  Rcpp::NumericMatrix Bd = Rcpp::as<Rcpp::NumericMatrix>(args_input["Bd"]);
  Rcpp::IntegerMatrix Bi = Rcpp::as<Rcpp::IntegerMatrix>(args_input["Bi"]);
  Rcpp::NumericVector B1d = Rcpp::as<Rcpp::NumericVector>(args_input["B1d"]);
  Rcpp::IntegerVector B1i = Rcpp::as<Rcpp::IntegerVector>(args_input["B1i"]);


  fmesh::Matrix<double> Bdd = Bd;
  fmesh::Matrix<int> Bdi(Bd);
  fmesh::Matrix<double> Bid(Bi);
  fmesh::Matrix<int> Bii(Bi);

  FMLOG_("Bdd: " << Bdd << std::endl);
  FMLOG_("Bdi: " << Bdi << std::endl);
  FMLOG_("Bid: " << Bid << std::endl);
  FMLOG_("Bii: " << Bii << std::endl);

  fmesh::Matrix1<double> Bdd1 = B1d;
  fmesh::Matrix1<double> Bdd1_ = Rcpp::NumericVector(Bd(Rcpp::_, 1));
  fmesh::Matrix1<double> Bdd1_0 = Bd;
  fmesh::Matrix1<int> Bdi1(B1d);
  fmesh::Matrix1<double> Bid1(B1i);
  fmesh::Matrix1<int> Bii1(B1i);

  FMLOG_("Bdd1: " << Bdd1 << std::endl);
  FMLOG_("Bdd1_: " << Bdd1_ << std::endl);
  FMLOG_("Bdd1_0: " << Bdd1_0 << std::endl);
  FMLOG_("Bdi1: " << Bdi1 << std::endl);
  FMLOG_("Bid1: " << Bid1 << std::endl);
  FMLOG_("Bii1: " << Bii1 << std::endl);

  fmesh::Matrix3<double> Bdd3 = Bd;
  fmesh::Matrix3<int> Bdi3(Bd);
  fmesh::Matrix3<double> Bid3(Bi);
  fmesh::Matrix3<int> Bii3(Bi);

  FMLOG_("Bdd3: " << Bdd3 << std::endl);
  FMLOG_("Bdi3: " << Bdi3 << std::endl);
  FMLOG_("Bid3: " << Bid3 << std::endl);
  FMLOG_("Bii3: " << Bii3 << std::endl);

  const EigenMSM<double> Ad(Rcpp::as<EigenMSM<double>>(args_input["Ad"]));

  fmesh::SparseMatrix<double> Ad_fm(Ad);

//  const EigenMSM<int> Ai(Rcpp::as<EigenMSM<int>>(args_input["Ai"]));

  //  bool is_msm = Rcpp::is<Eigen::SparseMatrix<double>>(args_input["a"]);

  matrices.attach("Ad_fm", &Ad_fm, false);
  matrices.output("Ad_fm");

  MatrixC mat2(args_input);

  FMLOG_(mat2.DI("tv"))

  Rcpp::List ret;
  ret["is_list"] = is_list;
  ret["is_numeric_matrix"] = is_numeric_matrix;
  ret["is_numeric_vector"] = is_numeric_vector;
  ret["is_integer_matrix"] = is_integer_matrix;
  ret["is_integer_vector"] = is_integer_vector;
//  ret["A"] = A;
  ret["Ad"] = Ad;
  ret["Ad_fm"] = Ad_fm.EigenSparseMatrix();
  ret["Ad_fm_auto"] = Ad_fm;
  ret["Ad_fm_ijx"] = Ad_fm.RcppList();
  ret["Bid3"] = Bid3;
  ret["Bdi3"] = Bdi3;
  ret["matrices"] = matrices;
  ret["mat2"] = mat2.output("-");
  return (ret);
}

#endif
