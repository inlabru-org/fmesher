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

#include "fmesher_helpers.h"
#include "fmesher.h"

using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
using std::endl;

using fmesh::Dart;
using fmesh::DartPair;
using fmesh::DartList;
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
using fmesh::Vector3;
using fmesh::constrMetaT;
using fmesh::constrT;
using fmesh::constrListT;
using fmesh::vertexListT;
using fmesh::TriangleLocator;

template<class T>
using EigenMSM = Eigen::MappedSparseMatrix<T>;
template<class T>
using EigenSM = Eigen::SparseMatrix<T>;

const bool useVT = true;
const bool useTTi = true;










#include "qtool.h"

//' Compute sparse matrix inverse
//'
//' @param AA A sparse matrix
// [[Rcpp::export]]
Rcpp::List C_qinv(SEXP AA)
{
  //Eigen::SparseMatrix<double> C_qinv(SEXP AA)
  const EigenMSM<double> A(Rcpp::as<EigenMSM<double> >(AA));

  QTool<double> Q;
  Q.Q(A);

  Rcpp::List ret;
  ret["Qinv"] = Q.S();
  return(ret);
  //  return Rcpp::List::create(Rcpp::Named("Q") = Q.S());
  //  return Q.S();
}




//' Triangulate
//'
//' @param args_input Input argument list
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_triangulate(Rcpp::List args_input)
{

  return Rcpp::List::create();
}

//' Test the matrix I/O system
//'
//' @param args_input Input argument list
// [[Rcpp::export]]
Rcpp::List C_matrixio_test(Rcpp::List args_input)
{
  //Eigen::SparseMatrix<double> C_qinv(SEXP AA)

  MatrixC matrices;

//  matrices.R_obj_input(input_matrix_list_from_R);
//  matrices.R_file_input(input_name_list_from_R);

  Rcpp::List ret;
  return(ret);
}



#endif
