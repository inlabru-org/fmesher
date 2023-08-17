#ifndef _FMESH_RCPP_FMESHER_
#define _FMESH_RCPP_FMESHER_ 1
#ifdef FMESHER_WITH_R

#include "Rcpp.h"
#ifdef FMESHER_WITH_EIGEN
#include "RcppEigen.h"
#endif

namespace fmesh {
template <class T> class Matrix;
template <class T> class Matrix1;
template <class T> class Matrix3;
template <class T> class SparseMatrix;
class MatrixC;

template <typename T>
struct Rcpp_traits;

template <>
struct Rcpp_traits<int> {
  using Vector = Rcpp::IntegerVector;
  using Matrix = Rcpp::IntegerMatrix;
};

template <>
struct Rcpp_traits<double> {
  using Vector = Rcpp::NumericVector;
  using Matrix = Rcpp::NumericMatrix;
};

/* Define default specialisation that can't match any real types. */
template <typename T>
struct Rcpp_traits {
private:
  struct VectorType;
  struct MatrixType;
public:
  using Vector = VectorType;
  using Matrix = MatrixType;
};

}

/* forward declarations */
namespace Rcpp {
/* support for wrap */

#define FM_DEFINE_WRAP(__thetype__)                     \
  template<> inline SEXP wrap(const fmesh::__thetype__& obj);

FM_DEFINE_WRAP(Matrix<double>);
FM_DEFINE_WRAP(Matrix<int>);
FM_DEFINE_WRAP(Matrix1<double>);
FM_DEFINE_WRAP(Matrix1<int>);
FM_DEFINE_WRAP(Matrix3<double>);
FM_DEFINE_WRAP(Matrix3<int>);
FM_DEFINE_WRAP(SparseMatrix<double>);
FM_DEFINE_WRAP(SparseMatrix<int>);
FM_DEFINE_WRAP(MatrixC);

// TODO:
/* support for as */
// template<typename T> class Exporter< fmesh::Matrix<T> >;

}

#endif
#endif
