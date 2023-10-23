#ifndef _FMESH_VECTOR_T_
#define _FMESH_VECTOR_T_ 1

#include "ioutils.h"
#include "vector.h"
#include <fstream>
#include <sstream>

namespace fmesh {

template <class T>
Matrix<T>::Matrix(size_t set_rows, size_t set_cols, const T *vals)
    : data_(NULL), rows_(0), cols_(0), cap_(0) {
  cols(set_cols);
  capacity(set_rows);
  rows_ = set_rows;
  if (vals) {
    std::memcpy(&data_[0], vals, sizeof(T) * rows_ * cols_);
  }
}

template <class T>
Matrix<T>::Matrix(const Matrix<T> &from)
    : data_(NULL), rows_(0), cols_(0), cap_(0) {
  clear();
  cols(from.cols_);
  capacity(from.cap_);
  rows_ = from.rows_;
  if (data_) {
    std::memcpy(&data_[0], &from.data_[0], sizeof(T) * rows_ * cols_);
  }
}
#ifdef FMESHER_WITH_R
template <class T>
Matrix<T>::Matrix(const typename Matrix<T>::RcppMatrix &from)
  : data_(NULL), rows_(0), cols_(0), cap_(0) {
  clear();
  cols(from.ncol());
  capacity(from.nrow());
  rows_ = from.nrow();
  if (data_) {
    for (size_t col_idx = 0; col_idx < cols_; col_idx++) {
      typename RcppMatrix::ConstColumn col = from(Rcpp::_, col_idx);
      size_t row_idx = 0;
      for (auto elem = col.cbegin(); elem != col.cend(); elem++) {
        (*this)(row_idx, col_idx) = (T)(*elem);
        row_idx++;
      }
    }
  }
}

/*
template <class T>
Matrix<T>::Matrix(const Rcpp::IntegerMatrix &from)
  : data_(NULL), rows_(0), cols_(0), cap_(0) {
  clear();
  cols(from.ncol());
  capacity(from.nrow());
  rows_ = from.nrow();
  if (data_) {
    for (size_t col_idx = 0; col_idx < cols_; col_idx++) {
      Rcpp::IntegerMatrix::ConstColumn col = from(Rcpp::_, col_idx);
      size_t row_idx = 0;
      for (auto elem = col.cbegin(); elem != col.cend(); elem++) {
        (*this)(row_idx, col_idx) = (T)(*elem);
        row_idx++;
      }
    }
  }
}
*/
template <class T>
Matrix<T>::Matrix(const typename Matrix<T>::RcppVector &from)
  : data_(NULL), rows_(0), cols_(0), cap_(0) {
  clear();
  cols(1);
  capacity(from.length());
  rows_ = from.length();
  if (data_) {
    size_t row_idx = 0;
    for (auto elem = from.cbegin(); elem != from.cend(); elem++) {
      (*this)(row_idx, 0) = (T)(*elem);
      row_idx++;
    }
  }
}

/*
template <class T>
Matrix<T>::Matrix(const Rcpp::IntegerVector &from)
  : data_(NULL), rows_(0), cols_(0), cap_(0) {
  clear();
  cols(1);
  capacity(from.length());
  rows_ = from.length();
  if (data_) {
    size_t row_idx = 0;
    for (auto elem = from.cbegin(); elem != from.cend(); elem++) {
      (*this)(row_idx, 0) = (T)(*elem);
      row_idx++;
    }
  }
}
 */
#endif


template <class T>
const Matrix<T> &Matrix<T>::operator=(const Matrix<T> &from) {
  clear();
  cols(from.cols_);
  capacity(from.cap_);
  rows_ = from.rows_;
  if (data_ && from.data_) {
    std::memcpy(&data_[0], &from.data_[0], sizeof(T) * rows_ * cols_);
  }
  return *this;
}

template <class T> bool Matrix<T>::capacity(size_t cap) {
  if (cap <= cap_) {
    return true;
  }
  size_t old_cap = cap_;
  if ((cap_ == 0) && (cap < capacity_step_size_))
    cap_ = cap;
  while (cap > cap_) {
    if (cap_ < capacity_step_size_)
      cap_ = capacity_step_size_;
    else {
      if (cap_ < capacity_doubling_limit_)
        cap_ *= 2;
      else
        cap_ += capacity_step_size_;
    }
  }

  std::unique_ptr<T[]> data_new_ = std::make_unique<T[]>(cap_ * cols_);

  if (data_ && data_new_) { /* Copy existing data: */
    std::memcpy(&data_new_[0], &data_[0], sizeof(T) * old_cap * cols_);
  }

  data_ = std::move(data_new_);
  zeros(old_cap, cap_);
  return true;
}

template <class T> bool Matrix<T>::append(const Matrix<T> &toappend) {
  if (cols_ > 0) {
    if (cols_ != toappend.cols_)
      return false;
  } else {
    cols(toappend.cols_);
  }
  if (!capacity(rows_ + toappend.rows_))
    return false;
  if (data_ && toappend.data_) {
    std::memcpy(&data_[0] + rows_ * cols_, &toappend.data_[0],
                sizeof(T) * toappend.rows_ * cols_);
    rows_ += toappend.rows_;
  }
  return true;
}

template <class T> Matrix<T> &Matrix<T>::rows(size_t set_rows) {
  if (set_rows > rows_) {
    capacity(set_rows);
  } else if (set_rows < rows_)
    zeros(set_rows, rows_);
  rows_ = set_rows;
  return *this;
}
template <class T> Matrix<T> &Matrix<T>::cols(size_t set_cols) {
  /* When set, cannot alter number of columns,
     unless we only need to expand a single row. */
  if ((cols_ > 0) && (!((rows_ <= 1) && cols_ <= set_cols))) {
    return *this;
  }
  if ((cols_ > 0) && (rows_ > 0)) { /* We already have some data, and
                                       need to carefully make sure the
                                       capacity is enough */
    /* Pre-data-size: cap_*cols_
       Post-data_size needed cap_*set_cols
       capacity is set as r*cols_
       Requirement: r*cols_ >= cap_*set_cols
       r = (cap_*set_cols)/cols_+1
    */
    capacity((cap_ * set_cols) / cols_ + 1);
    cols_ = set_cols;
    cap_ = rows_; /* This makes sure we dont' overestimate the true cap_. */
  } else
    cols_ = set_cols;
  return *this;
}



#ifdef FMESHER_WITH_R
template <class T>
Matrix1<T>::Matrix1(const typename Matrix<T>::RcppMatrix &from) : Matrix<T>(from) {
  // Check number of columns:
  if (1 != from.ncol()) {
    Rcpp::stop("NumericMatrix->Matrix1 conversion: column number mismatch.");
  }
}

template <class T>
Matrix1<T>::Matrix1(const typename Matrix<T>::RcppVector &from) : Matrix<T>(from) {}

template <class T>
Matrix3<T>::Matrix3(const typename Matrix<T>::RcppMatrix &from) : Matrix<T>(from) {
  // Check number of columns:
  if (3 != from.ncol()) {
    Rcpp::stop("NumericMatrix->Matrix3 conversion: column number mismatch.");
  }
}

#endif




// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R

template <class T>
bool Matrix<T>::save(std::string filename, IOMatrixtype matrixt,
                     bool binary) const {
  std::ofstream O;
  if (binary)
    O.open(filename.c_str(), std::ios::out | std::ios::binary);
  else
    O.open(filename.c_str(), std::ios::out);
  if (!O.is_open())
    return false;
  IOHelperM<T> ioh;
  ioh.cD(this).binary(binary).matrixtype(matrixt);
  ioh.OH(O).OD(O);
  O.close();
  return true;
}

template <class T> bool Matrix<T>::load(std::string filename, bool binary) {
  std::ifstream I;
  if (binary)
    I.open(filename.c_str(), std::ios::in | std::ios::binary);
  else
    I.open(filename.c_str(), std::ios::in);
  if (!I.is_open())
    return false;
  IOHelperM<T> ioh;
  ioh.D(this);
  ioh.binary(binary).IH(I).ID(I);
  I.close();
  return true;
}

template <class T>
bool Matrix<T>::save_ascii_2009(std::string filename,
                                IOMatrixtype matrixt) const {
  std::ofstream O;
  O.open(filename.c_str(), std::ios::out);
  if (!O.is_open())
    return false;
  IOHelperM<T> ioh;
  ioh.cD(this).ascii().matrixtype(matrixt);
  ioh.OD_2009(O);
  O.close();
  return true;
}

template <class T> bool Matrix<T>::load_ascii_2009(std::string filename) {
  std::ifstream I;
  I.open(filename.c_str(), std::ios::in);
  if (!I.is_open())
    return false;
  load_ascii_2009(I);
  I.close();
  return true;
}

template <class T> void Matrix<T>::load_ascii_2009(std::istream &input) {
  (*this).clear();
  int r = 0;
  int cols = 0;
  while (!input.eof()) {
    std::string line;
    getline(input, line);
    std::stringstream ss(line);
    int c = 0;
    if (line.length() > 0)
      while (!ss.eof()) {
        ss >> (*this)(r, c);
        c++;
      }
    if (r == 0) {
      cols = c;
    }
    if (c < cols) {
      (*this).rows(r);
      r--;
    }
    r++;
  }
}

template <class T>
bool SparseMatrix<T>::save(std::string filename, IOMatrixtype matrixt,
                           bool binary) const {
  std::ofstream O;
  if (binary)
    O.open(filename.c_str(), std::ios::out | std::ios::binary);
  else
    O.open(filename.c_str(), std::ios::out);
  if (!O.is_open())
    return false;
  IOHelperSM<T> ioh;
  ioh.cD(this).binary(binary).matrixtype(matrixt);
  ioh.OH(O).OD(O);
  O.close();
  return true;
}

template <class T>
bool SparseMatrix<T>::load(std::string filename, bool binary) {
  std::ifstream I;
  if (binary)
    I.open(filename.c_str(), std::ios::in | std::ios::binary);
  else
    I.open(filename.c_str(), std::ios::in);
  if (!I.is_open())
    return false;
  IOHelperSM<T> ioh;
  ioh.D(this);
  ioh.binary(binary).IH(I).ID(I);
  I.close();
  return true;
}
#endif // not FMESHER_WITH_R

#ifdef FMESHER_WITH_R

// Export to R
#ifdef FMESHER_WITH_EIGEN
template <class T> Eigen::SparseMatrix<T> SparseMatrix<T>::EigenSparseMatrix(IOMatrixtype matrixt) const {
  typedef Eigen::Triplet<T> Trip;
  std::vector<Trip> tripletList;
  tripletList.reserve(nnz());
  tolist(tripletList, matrixt);

  Eigen::SparseMatrix<T> m(rows(), cols());
  m.setFromTriplets(tripletList.begin(), tripletList.end());

  return m;
}
#endif

// Export to R as list(i,j,x,dims)
template <class T> SEXP SparseMatrix<T>::fmesher_sparse(IOMatrixtype matrixt) const {
  std::vector<int> i;
  std::vector<int> j;
  std::vector<T> x;
  std::vector<int> dims;
  to_ijx(i, j, x, dims, matrixt);

  Rcpp::List output;
  output["i"] = i;
  output["j"] = j;
  output["x"] = x;
  output["dims"] = dims;
  output.attr("class") = "fmesher_sparse";

  return output;
}

template <class T> SEXP SparseMatrix<T>::unpackedMatrix(IOMatrixtype matrixt) const {
  Rcpp::List output = fmesher_sparse(matrixt);
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("fmesher");
  Rcpp::Function fun = pkg["fm_as_unpackedMatrix"];
  Rcpp::S4 obj = fun(output);

  return obj;
}

template <class T> SEXP SparseMatrix<T>::dgCMatrix(IOMatrixtype matrixt) const {
  Rcpp::List output = fmesher_sparse(matrixt);
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("fmesher");
  Rcpp::Function fun = pkg["fm_as_dgCMatrix"];
  Rcpp::S4 obj = fun(output);

  return obj;
}

template <class T> SEXP SparseMatrix<T>::dgTMatrix(IOMatrixtype matrixt) const {
  Rcpp::List output = fmesher_sparse(matrixt);
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("fmesher");
  Rcpp::Function fun = pkg["fm_as_dgTMatrix"];
  Rcpp::S4 obj = fun(output);

  return obj;
}


template <class T>
void SparseMatrix<T>::fromRcpp(SEXP from) {
  if (Rcpp::is<Rcpp::List>(from) &&
      Rcpp::as<Rcpp::List>(from).inherits("fmesher_sparse")) {
    Rcpp::List from_list = Rcpp::as<Rcpp::List>(from);
    Rcpp::IntegerVector Tr = Rcpp::as<Rcpp::IntegerVector>(from_list["i"]);
    Rcpp::IntegerVector Tc = Rcpp::as<Rcpp::IntegerVector>(from_list["j"]);
    Rcpp::NumericVector Tv = Rcpp::as<Rcpp::NumericVector>(from_list["x"]);
    Rcpp::IntegerVector dims = Rcpp::as<Rcpp::IntegerVector>(from_list["dims"]);

    fromlist(Tr, Tc, Tv, dims, IOMatrixtype::General);
  } else if (Rcpp::is<Rcpp::S4>(from)) {
    Rcpp::S4 obj = (SEXP)from;
    if (obj.is("Matrix")) {
      if (!obj.is("dgTMatrix")) {
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("fmesher");
        Rcpp::Function fun = pkg["fm_as_dgTMatrix"];
        obj = fun(from);
      }

      Rcpp::IntegerVector Tr = Rcpp::as<Rcpp::IntegerVector>(obj.slot("i"));
      Rcpp::IntegerVector Tc = Rcpp::as<Rcpp::IntegerVector>(obj.slot("j"));
      Rcpp::NumericVector Tv = Rcpp::as<Rcpp::NumericVector>(obj.slot("x"));
      Rcpp::IntegerVector dims = Rcpp::as<Rcpp::IntegerVector>(obj.slot("Dim"));

      fromlist(Tr, Tc, Tv, dims, IOMatrixtype::General);
    } else {
      Rcpp::warning("Unsupported SparseMatrix<T>(Rcpp::S4) class.");
    }
  } else {
    Rcpp::warning("Unsupported SparseMatrix<T>(Rcpp) class.");
  }
}

template <class T>
SparseMatrix<T>::SparseMatrix(SEXP from)
  : cols_(0), data_() {
  fromRcpp(from);
}

#ifdef FMESHER_WITH_EIGEN

template<class T> using EigenMSM = Eigen::Map<Eigen::SparseMatrix<T>>;

template <class T>
void SparseMatrix<T>::fromEigen(const EigenMSM<T> &from) {
  rows(from.rows());
  for (int k=0; k < from.outerSize(); ++k)
    for (typename EigenMSM<T>::InnerIterator it(from, k); it; ++it)
    {
      operator()(it.row(), it.col(), it.value());
    }
}

template <class T>
SparseMatrix<T>::SparseMatrix(const EigenMSM<T> &from)
  : cols_(from.cols()), data_() {
  fromEigen(from);
}

#endif

#endif


// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R

template <class T>
bool SparseMatrix<T>::save_ascii_2009(std::string filename,
                                      IOMatrixtype matrixt) const {
  std::ofstream O;
  O.open(filename.c_str(), std::ios::out);
  if (!O.is_open())
    return false;
  IOHelperSM<T> ioh;
  ioh.cD(this).ascii().matrixtype(matrixt);
  ioh.OD_2009(O);
  O.close();
  return true;
}

#endif // not FMESHER_WITH_R

template <class T> SparseMatrix<T> diag(const Matrix<T> &M1) {
  SparseMatrix<T> SM(M1.rows(), M1.rows());
  for (size_t i = 0; i < M1.rows(); i++) {
    SM(i, i) = M1[i][0];
  }
  return SM;
}

template <class T> Matrix<T> diag(const SparseMatrix<T> &M1) {
  Matrix<T> M(M1.rows(), 1);
  for (size_t i = 0; ((i < M1.rows()) && (i < M1.cols())); i++) {
    M(i, 0) = M1[i][i];
  }
  return M;
}

template <class T>
SparseMatrix<T> operator*(const SparseMatrix<T> &M1,
                          const SparseMatrix<T> &M2) {
  SparseMatrix<T> M;
  size_t M1rows = M1.rows();
  size_t M2rows = M2.rows();
  M.cols(M2.cols()).rows(M1rows);
  for (size_t i = 0; i < M1rows; i++) {
    SparseMatrixRow<T> &Mi = M(i);
    const SparseMatrixRow<T> &M1i = M1[i];
    if (M1i.size() > 0) {
      for (typename SparseMatrixRow<T>::ColCIter M1k = M1i.begin();
           (M1k != M1i.end()) && ((size_t)M1k->first < M2rows); M1k++) {
        int k = M1k->first;
        const T &M1ik = M1i[k];
        const SparseMatrixRow<T> &M2k = M2[k];
        for (typename SparseMatrixRow<T>::ColCIter M2j = M2k.begin();
             (M2j != M2k.end()); M2j++) {
          Mi(M2j->first) += (M1ik * M2j->second);
        }
      }
    }
  }
  return M;
}

template <class T>
SparseMatrix<T> operator-(const SparseMatrix<T> &M1,
                          const SparseMatrix<T> &M2) {
  SparseMatrix<T> M(M1);
  for (size_t r = 0; (r < M1.rows()) && (r < M2.rows()); r++) {
    SparseMatrixRow<T> &Mr = M(r);
    const SparseMatrixRow<T> &M2r = M2[r];
    for (typename SparseMatrixRow<T>::ColCIter c = M2r.begin();
         (c != M2r.end()) && ((size_t)c->first < M1.cols()); c++) {
      Mr(c->first) -= c->second;
    }
  }
  return M;
}

template <class T>
SparseMatrix<T> inverse(const SparseMatrix<T> &M1, bool diagonal) {
  SparseMatrix<T> M;
  M.cols(M1.cols()).rows(M1.rows());
  if (!diagonal) {
    /* NOT IMPLEMENTED */
    NOT_IMPLEMENTED;
    return M;
  }
  for (size_t r = 0; (r < M1.rows()) && (r < M1.cols()); r++) {
    const T &val = M1[r][r];
    if (!(val == T()))
      M(r, r) = 1 / val;
  }
  return M;
}

template <class T> const T fmesh::Matrix<T>::zero_ = T();
template <class T> const T fmesh::SparseMatrixRow<T>::zero_ = T();
template <class T> const T fmesh::SparseMatrix<T>::zero_ = T();

} /* namespace fmesh */

#ifdef FMESHER_WITH_R
namespace Rcpp {

#define __FM_MATRIX_WRAP__(InType, T, COLS)                    \
template<>                                                     \
inline SEXP wrap(const fmesh::InType<T>& obj) {                \
  fmesh::Rcpp_traits<T>::Matrix res(obj.rows(), COLS);         \
  for (size_t r = 0; r < obj.rows(); r++) {                    \
    for (size_t c = 0; c < COLS; c++) {                        \
      res(r, c) = obj[r][c];                                   \
    }                                                          \
  }                                                            \
  return res;                                                  \
}
#define __FM_VECTOR_WRAP__(InType, T)                          \
template<>                                                     \
inline SEXP wrap(const fmesh::InType<T>& obj) {                \
  fmesh::Rcpp_traits<T>::Vector res(obj.rows());               \
  for (size_t r = 0; r < obj.rows(); r++) {                    \
    res[r] = obj[r];                                           \
  }                                                            \
  return res;                                                  \
}

__FM_MATRIX_WRAP__(Matrix, double, obj.cols())
__FM_MATRIX_WRAP__(Matrix, int, obj.cols())
__FM_MATRIX_WRAP__(Matrix3, double, 3)
__FM_MATRIX_WRAP__(Matrix3, int, 3)
__FM_VECTOR_WRAP__(Matrix1, double)
__FM_VECTOR_WRAP__(Matrix1, int)

  template<>
  inline SEXP wrap(const fmesh::SparseMatrix<double>& obj) {
    return Rcpp::wrap(obj.dgTMatrix(fmesh::IOMatrixtype::General));
  }
  template<>
  inline SEXP wrap(const fmesh::SparseMatrix<int>& obj) {
    // No Sparse matrix storage for integers, so return ijx triples instead.
    return Rcpp::wrap(obj.fmesher_sparse(fmesh::IOMatrixtype::General));
  }

} // Namespace Rcpp
#endif

#endif
