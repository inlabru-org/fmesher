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

#include "Rcpp.h"
#include "RcppEigen.h"

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

const bool useVT = true;
const bool useTTi = true;

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

//' Triangulate
//'
//' @param args_input Input argument list
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_triangulate(Rcpp::List args_input) {

  return Rcpp::List::create();
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
//  const EigenMSM<int> Ai(Rcpp::as<EigenMSM<int>>(args_input["Ai"]));

  //  bool is_msm = Rcpp::is<Eigen::SparseMatrix<double>>(args_input["a"]);

  //  matrices.R_obj_input(input_matrix_list_from_R);
  //  matrices.R_file_input(input_name_list_from_R);

  Rcpp::List ret;
  ret["is_list"] = is_list;
  ret["is_numeric_matrix"] = is_numeric_matrix;
  ret["is_numeric_vector"] = is_numeric_vector;
  ret["is_integer_matrix"] = is_integer_matrix;
  ret["is_integer_vector"] = is_integer_vector;
//  ret["A"] = A;
  ret["Ad"] = Ad;
  //  ret["is_msm"] = is_msm;
  return (ret);
}

#endif
