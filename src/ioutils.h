#ifndef _FMESH_IOUTILS_
#define _FMESH_IOUTILS_ 1

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <variant>
#include <vector>

#include "fmesher_debuglog.h"
#include "vector.h"

#define IOHEADER_VERSION 0
#define BINARY_DEFAULT false

namespace fmesh {

class MatrixC;
class IOHeader;
template <class T> class IOHelper;
template <class T> class IOHelperM;
template <class T> class IOHelperSM;
class IOHelperC;

/*! dense/sparse/map */
enum class IODatatype : int {
  invalid = -1,
    dense = 0,
    sparse = 1,
    collection = 2
};
/*! int/double */
enum IOValuetype { IOValuetype_int = 0, IOValuetype_double = 1 };
/*! rowmajor/colmajor */
enum IOStoragetype { IOStoragetype_rowmajor = 0, IOStoragetype_colmajor = 1 };


// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R

/*! Header for input and output file formats. */
class IOHeader {
public:
  int version;     /*!< Format version */
  int elems;       /*!< The number of data units

                     For dense matrices, the total number of elements.
                     For sparse matrices, the number of elements contained in
                     the file.
                    */
  int rows;        /*!< The number of data rows. */
  int cols;        /*!< The number of data columns. */
  IODatatype datatype;    /*!< The IODatatype. */
  int valuetype;   /*!< The IOValuetype. */
  int matrixtype;  /*!< The IOMatrixtype. */
  int storagetype; /*!< The IOStoragetype. */

  /* Sets defaults, and the valuetype matching T: */
  template <class T> IOHeader &def(const T &ref);
  IOHeader &def(const int &ref);
  IOHeader &def(const double &ref);
  IOHeader &def(const MatrixC &ref);
  IOHeader &def();
  /* Default values: */
  template <class T>
  IOHeader &dense(const Matrix<T> &M,
                  IOMatrixtype matrixt = IOMatrixtype_general);
  template <class T>
  IOHeader &sparse(const SparseMatrix<T> &M,
                   IOMatrixtype matrixt = IOMatrixtype_general);
  IOHeader &collection(const MatrixC &C);

  /* Constructor, that sets the valuetype matching T: */
  template <class T> IOHeader(const T &ref);
  IOHeader();
};

std::ostream &operator<<(std::ostream &output, const IOHeader &h);
std::istream &operator>>(std::istream &output, IOHeader &h);

/*! Base helper for input and output. */
template <class T> class IOHelper {
public:
  IOHeader h_;
  bool binary_;

public:
  /* Constructors: */
  IOHelper() : h_(T()), binary_(BINARY_DEFAULT){};
  IOHelper(const IOHeader &h) : h_(h), binary_(BINARY_DEFAULT){};

  bool binaryformat() const { return binary_; };
  IOMatrixtype matrixtype() const { return (IOMatrixtype)h_.matrixtype; };
  IOStoragetype storagetype() const { return (IOStoragetype)h_.storagetype; };

  IOHelper<T> &ascii(bool set_ascii = true) { return binary(!set_ascii); };
  IOHelper<T> &binary(bool set_binary = true) {
    binary_ = set_binary;
    return *this;
  };
  IOHelper<T> &storagetype(IOStoragetype set_storage) {
    h_.storagetype = set_storage;
    return *this;
  };
  IOHelper<T> &rowmajor(bool set_rowmajor = true) {
    h_.storagetype =
        (set_rowmajor ? IOStoragetype_rowmajor : IOStoragetype_rowmajor);
    return *this;
  };
  IOHelper<T> &colmajor(bool set_colmajor = true) {
    h_.storagetype =
        (set_colmajor ? IOStoragetype_colmajor : IOStoragetype_colmajor);
    return *this;
  };

  /* Output/Input: */
  IOHelper<T> &OH(std::ostream &output);
  IOHelper<T> &IH(std::istream &input);
  IOHelper<T> &IH(const IOHeader &h);
};

/*! Helper for Matrix input and output. */
template <class T> class IOHelperM : public IOHelper<T> {
public:
  const Matrix<T> *cM_;
  Matrix<T> *M_;

public:
  /* Constructors: */
  IOHelperM() : IOHelper<T>(), cM_(NULL), M_(NULL) {};
  IOHelperM(const IOHeader &h) : IOHelper<T>(h), cM_(NULL), M_(NULL) {};
  IOHelperM<T> &cD(const Matrix<T> *M) {
    cM_ = M;
    M_ = NULL;
    IOHelper<T>::h_.dense(*M);
    return *this;
  };
  IOHelperM<T> &D(Matrix<T> *M) {
    cM_ = M;
    M_ = M;
    IOHelper<T>::h_.dense(*M);
    return *this;
  };

  IOHelperM<T> &matrixtype(IOMatrixtype matrixt) {
    IOHelper<T>::h_.dense(*cM_, matrixt);
    return *this;
  };

  /* Output/Input: */
  IOHelperM<T> &OD(std::ostream &output);
  IOHelperM<T> &ID(std::istream &input);

  /* Backwards compatibility: */
  IOHelperM<T> &OH_2009(std::ostream &output);
  IOHelperM<T> &OD_2009(std::ostream &output);

  /* Overloaded from IOHelper: */
  IOHelperM<T> &ascii(bool set_ascii = true) {
    IOHelper<T>::ascii(set_ascii);
    return *this;
  };
  IOHelperM<T> &binary(bool set_binary = true) {
    IOHelper<T>::binary(set_binary);
    return *this;
  };
  IOHelperM<T> &general() { return matrixtype(IOMatrixtype_general); };
  IOHelperM<T> &symmetric() { return matrixtype(IOMatrixtype_symmetric); };
  IOHelperM<T> &diagonal() { return matrixtype(IOMatrixtype_diagonal); };
  IOHelperM<T> &storagetype(IOStoragetype set_storage) {
    IOHelper<T>::storagetype(set_storage);
    return *this;
  };
  IOHelperM<T> &rowmajor(bool set_rowmajor = true) {
    IOHelper<T>::rowmajor(set_rowmajor);
    return *this;
  };
  IOHelperM<T> &colmajor(bool set_colmajor = true) {
    IOHelper<T>::colmajor(set_colmajor);
    return *this;
  };
  IOHelperM<T> &OH(std::ostream &output) {
    IOHelper<T>::OH(output);
    return *this;
  };
  IOHelperM<T> &IH(std::istream &input) {
    IOHelper<T>::IH(input);
    return *this;
  };
  IOHelperM<T> &IH(const IOHeader &h) {
    IOHelper<T>::IH(h);
    return *this;
  };
};

/*! Helper for SparseMatrix input and output. */
template <class T> class IOHelperSM : public IOHelper<T> {
public:
  const SparseMatrix<T> *cM_;
  SparseMatrix<T> *M_;

public:
  /* Constructors: */
  IOHelperSM() : IOHelper<T>(), cM_(NULL), M_(NULL){};
  IOHelperSM(const IOHeader &h) : IOHelper<T>(h), cM_(NULL), M_(NULL){};
  IOHelperSM<T> &cD(const SparseMatrix<T> *M) {
    cM_ = M;
    M_ = NULL;
    IOHelper<T>::h_.sparse(*M);
    IOHelper<T>::colmajor();
    return *this;
  };
  IOHelperSM<T> &D(SparseMatrix<T> *M) {
    cM_ = M;
    M_ = M;
    IOHelper<T>::h_.sparse(*M);
    IOHelper<T>::colmajor();
    return *this;
  };

  IOHelperSM<T> &matrixtype(IOMatrixtype matrixt) {
    IOHelper<T>::h_.sparse(*cM_, matrixt);
    return *this;
  };

  /* Output/Input: */
  IOHelperSM<T> &OD(std::ostream &output);
  IOHelperSM<T> &ID(std::istream &input);
  /* Backwards compatibility: */
  IOHelperSM<T> &OH_2009(std::ostream &output);
  IOHelperSM<T> &OD_2009(std::ostream &output);

  /* Overloaded from IOHelper: */
  IOHelperSM<T> &ascii(bool set_ascii = true) {
    binary(!set_ascii);
    return *this;
  };
  IOHelperSM<T> &binary(bool set_binary = true) {
    IOHelper<T>::binary(set_binary);
    if (set_binary)
      colmajor();
    else
      rowmajor();
    return *this;
  };
  IOHelperSM<T> &general() { return matrixtype(IOMatrixtype_general); };
  IOHelperSM<T> &symmetric() { return matrixtype(IOMatrixtype_symmetric); };
  IOHelperSM<T> &diagonal() { return matrixtype(IOMatrixtype_diagonal); };
  IOHelperSM<T> &storagetype(IOStoragetype set_storage) {
    IOHelper<T>::storagetype(set_storage);
    return *this;
  };
  IOHelperSM<T> &rowmajor(bool set_rowmajor = true) {
    IOHelper<T>::rowmajor(set_rowmajor);
    return *this;
  };
  IOHelperSM<T> &colmajor(bool set_colmajor = true) {
    IOHelper<T>::colmajor(set_colmajor);
    return *this;
  };
  IOHelperSM<T> &OH(std::ostream &output) {
    IOHelper<T>::OH(output);
    return *this;
  };
  IOHelperSM<T> &IH(std::istream &input) {
    IOHelper<T>::IH(input);
    return *this;
  };
  IOHelperSM<T> &IH(const IOHeader &h) {
    IOHelper<T>::IH(h);
    return *this;
  };
};

/*! Helper for MatrixC input and output. */
class IOHelperC : public IOHelper<int> {
public:
  typedef std::vector<std::string> listT;
  const MatrixC *cM_;
  MatrixC *M_;
  listT list_;

public:
  /* Constructors: */
  IOHelperC() : IOHelper<int>(IOHeader()), cM_(NULL), M_(NULL){};
  IOHelperC(const IOHeader &h) : IOHelper<int>(h), cM_(NULL), M_(NULL){};
  IOHelperC &cD(const MatrixC *M) {
    cM_ = M;
    M_ = NULL;
    IOHelper<int>::h_.collection(*M);
    return *this;
  };
  IOHelperC &D(MatrixC *M) {
    cM_ = M;
    M_ = M;
    IOHelper<int>::h_.collection(*M);
    return *this;
  };

  /* Output/Input: */
  IOHelperC &OL(std::ostream &output);
  IOHelperC &IL(std::istream &input);
  IOHelperC &OD(std::ostream &output);
  IOHelperC &ID(std::istream &input);

  /* Overloaded from IOHelper: */
  IOHelperC &ascii(bool set_ascii = true) {
    IOHelper<int>::ascii(set_ascii);
    return *this;
  };
  IOHelperC &binary(bool set_binary = true) {
    IOHelper<int>::binary(set_binary);
    return *this;
  };
  IOHelperC &OH(std::ostream &output) {
    IOHelper<int>::OH(output);
    return *this;
  };
  IOHelperC &IH(std::istream &input) {
    IOHelper<int>::IH(input);
    return *this;
  };
  IOHelperC &IH(const IOHeader &h) {
    IOHelper<int>::IH(h);
    return *this;
  };
};

#endif // not FMESHER_WITH_R


class MCCInfo {
public:
  bool loaded;
  bool active;
  IODatatype datatype;
  IOValuetype valuetype;
  IOMatrixtype matrixtype;
  bool owner;

  MCCInfo()
      : loaded(false), active(false), datatype(IODatatype::dense),
        valuetype(IOValuetype_int), matrixtype(IOMatrixtype_general),
        owner(false){};
  MCCInfo(bool load, bool act, IODatatype data, IOValuetype value,
          IOMatrixtype matrixt, bool isowner)
      : loaded(load), active(act), datatype(data), valuetype(value),
        matrixtype(matrixt), owner(isowner){};
};

typedef
std::variant<
  std::monostate,
  std::unique_ptr<Matrix<int>>,
  std::unique_ptr<Matrix<double>>,
  std::unique_ptr<SparseMatrix<int>>,
  std::unique_ptr<SparseMatrix<double>>,
  Matrix<int>*,
  Matrix<double>*,
  SparseMatrix<int>*,
  SparseMatrix<double>*
> MatrixVariantPtr;

class MCC {
  friend class MatrixC;

public:
  MCCInfo info;

protected:
  MatrixVariantPtr matrix_;

public:
  MCC()
      : info(false, false, IODatatype::dense, IOValuetype_int,
             IOMatrixtype_general, false),
        matrix_() {};
  MCC(IODatatype data, IOValuetype value, IOMatrixtype matrixt)
    : info(true, false, data, value, matrixt, true),
      matrix_() {
    if (info.datatype == IODatatype::dense) {
      if (info.valuetype == IOValuetype_int) {
        matrix_ = std::make_unique<Matrix<int>>();
      } else {
        matrix_ = std::make_unique<Matrix<double>>();
      }
    } else {
      if (info.valuetype == IOValuetype_int) {
        matrix_ = std::make_unique<SparseMatrix<int>>();
      } else {
        matrix_ = std::make_unique<SparseMatrix<double>>();
      }
    }
    if (std::holds_alternative<std::unique_ptr<Matrix<int>>>(matrix_)) {
      info.datatype = IODatatype::dense;
      info.valuetype = IOValuetype_int;
    } else if (std::holds_alternative<std::unique_ptr<Matrix<double>>>(matrix_)) {
      info.datatype = IODatatype::dense;
      info.valuetype = IOValuetype_double;
    } else if (std::holds_alternative<std::unique_ptr<SparseMatrix<int>>>(matrix_)) {
      info.datatype = IODatatype::sparse;
      info.valuetype = IOValuetype_int;
    } else if (std::holds_alternative<std::unique_ptr<SparseMatrix<double>>>(matrix_)) {
      info.datatype = IODatatype::sparse;
      info.valuetype = IOValuetype_double;
    }
  };
  template <class MatrixType>
  MCC(IODatatype data, IOValuetype value, IOMatrixtype matrixt,
      MatrixType * M,
      bool isowner = true)
    : info(true, false, data, value, matrixt, isowner), matrix_() {
    matrix_ = M;
    if (std::holds_alternative<Matrix<int>*>(matrix_)) {
      info.datatype = IODatatype::dense;
      info.valuetype = IOValuetype_int;
    } else if (std::holds_alternative<Matrix<double>*>(matrix_)) {
      info.datatype = IODatatype::dense;
      info.valuetype = IOValuetype_double;
    } else if (std::holds_alternative<SparseMatrix<int>*>(matrix_)) {
      info.datatype = IODatatype::sparse;
      info.valuetype = IOValuetype_int;
    } else if (std::holds_alternative<SparseMatrix<double>*>(matrix_)) {
      info.datatype = IODatatype::sparse;
      info.valuetype = IOValuetype_double;
    }
    if (isowner) {
      matrix_ = std::unique_ptr<MatrixType>(M);
    }
  };
  ~MCC() {
  };

  auto &DI() {
    if (std::holds_alternative<std::unique_ptr<Matrix<int>>>(matrix_)) {
      return *std::get<std::unique_ptr<Matrix<int>>>(matrix_);
    } else {
      return *std::get<Matrix<int>*>(matrix_);
    }
  }
  auto &DD() {
    if (std::holds_alternative<std::unique_ptr<Matrix<double>>>(matrix_)) {
      return *std::get<std::unique_ptr<Matrix<double>>>(matrix_);
    } else {
      return *std::get<Matrix<double>*>(matrix_);
    }
  }
  auto &SI() {
    if (std::holds_alternative<std::unique_ptr<SparseMatrix<int>>>(matrix_)) {
      return *std::get<std::unique_ptr<SparseMatrix<int>>>(matrix_);
    } else {
      return *std::get<SparseMatrix<int>*>(matrix_);
    }
  }
  auto &SD() {
    if (std::holds_alternative<std::unique_ptr<SparseMatrix<double>>>(matrix_)) {
      return *std::get<std::unique_ptr<SparseMatrix<double>>>(matrix_);
    } else {
      return *std::get<SparseMatrix<double>*>(matrix_);
    }
  }
  const auto &DI() const {
    if (std::holds_alternative<std::unique_ptr<Matrix<int>>>(matrix_)) {
      return *std::get<std::unique_ptr<Matrix<int>>>(matrix_);
    } else {
      return *std::get<Matrix<int>*>(matrix_);
    }
  }
  const auto &DD() const {
    if (std::holds_alternative<std::unique_ptr<Matrix<double>>>(matrix_)) {
      return *std::get<std::unique_ptr<Matrix<double>>>(matrix_);
    } else {
      return *std::get<Matrix<double>*>(matrix_);
    }
  }
  const auto &SI() const {
    if (std::holds_alternative<std::unique_ptr<SparseMatrix<int>>>(matrix_)) {
      return *std::get<std::unique_ptr<SparseMatrix<int>>>(matrix_);
    } else {
      return *std::get<SparseMatrix<int>*>(matrix_);
    }
  }
  const auto &SD() const {
    if (std::holds_alternative<std::unique_ptr<SparseMatrix<double>>>(matrix_)) {
      return *std::get<std::unique_ptr<SparseMatrix<double>>>(matrix_);
    } else {
      return *std::get<SparseMatrix<double>*>(matrix_);
    }
  }
};

class MatrixC {
  // No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
  friend class IOHelperC;
#endif
  typedef std::pair<std::string, MCC *> collPairT;
  typedef std::map<std::string, MCC *> collT;
  typedef std::set<std::string> outputT;
  typedef std::map<std::string, std::string> sourceT;

  collT coll_; /* name --> matrixdata */
  bool output_all_;
  outputT output_; /* names */
  bool bin_in_;
  bool bin_out_;
  sourceT source_; /* name --> filename */
  std::string input_prefix_;
  std::string output_prefix_;
  std::string output_file_;

public:
  MatrixC()
      : output_all_(false), bin_in_(true), bin_out_(true), input_prefix_("-"),
        output_prefix_("-"), output_file_(""){};
#ifdef FMESHER_WITH_R
  MatrixC(SEXP from);
  void attach(SEXP from); // Rccp::List of matrices
  void attach(std::string name, SEXP from);
#endif
  ~MatrixC() {
    for (auto& colli : coll_) {
      delete colli.second;
    }
  };

  int output_size() const { return output_.size(); }
  MatrixC &dont_output(std::string name);
  MatrixC &output(std::string name);

  void io(bool bin_in, bool bin_out);
  void input_prefix(std::string prefix);
  void output_prefix(std::string prefix);

// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
  void input_file(std::string filename);
  void output_file(std::string filename);

  template <class T> void input_raw_M(std::istream &input, Matrix<T> &M) const;

  void input_raw(std::string name, std::string specification,
                 std::string filename);
  void save();
#endif
#ifdef FMESHER_WITH_R
  SEXP Rcpp_wrap() const;
#endif

  // No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
  void load_file(std::string filename, bool only_list = false);
#endif

  /*! Activate all loaded matrices */
  void activate();
  /*! Activate if loaded */
  bool activate(std::string name);

// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
  /*! Load and activate */
  MCCInfo load(std::string name);
#endif

  /*! Add and activate */
  template <class T>
  Matrix<T> &attach(
      std::string name,
      Matrix<T> *M,
      bool transfer_ownership = true,
      IOMatrixtype matrixt = IOMatrixtype_general);
  template <class T>
  SparseMatrix<T> &attach(
      std::string name,
      SparseMatrix<T> *M,
      bool transfer_ownership = true,
      IOMatrixtype matrixt = IOMatrixtype_general);

  MatrixC &free(std::string name);

  Matrix<int> &DI(std::string name);
  Matrix<double> &DD(std::string name);
  SparseMatrix<int> &SI(std::string name);
  SparseMatrix<double> &SD(std::string name);

  void matrixtype(std::string name, IOMatrixtype matrixt);

  MCCInfo info(std::string name) const;
};

} /* namespace fmesh */

#include "ioutils_t.h"

#endif
