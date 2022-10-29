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



#include "qtool.h"

//' Compute sparse matrix inverse
//'
//' @param AA A sparse matrix
// [[Rcpp::export]]
Rcpp::List C_qinv(SEXP AA)
{
  //Eigen::SparseMatrix<double> C_qinv(SEXP AA)
  using Eigen::MappedSparseMatrix;
  using Eigen::SparseMatrix;
  const MappedSparseMatrix<double> A(Rcpp::as<MappedSparseMatrix<double> >(AA));

  QTool<double> Q;
  //  Q.Q(Rcpp::as<MappedSparseMatrix<double> >(AA));
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
  using Eigen::MappedSparseMatrix;
  using Eigen::SparseMatrix;

  Rcpp::List ret;
  return(ret);
}
