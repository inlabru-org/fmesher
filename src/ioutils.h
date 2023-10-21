#ifndef _FMESH_IOUTILS_
#define _FMESH_IOUTILS_ 1

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
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
  Invalid = -1,
    Dense = 0,
    Sparse = 1,
    Collection = 2
};
/*! int/double */
enum class IOValuetype : int {
  Invalid = -1,
    Int = 0,
    Double = 1
};
/*! rowmajor/colmajor */
enum class IOStoragetype : int {
  Invalid = -1,
    Rowmajor = 0,
    Colmajor = 1
};


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
  IOValuetype valuetype;   /*!< The IOValuetype. */
  IOMatrixtype matrixtype;  /*!< The IOMatrixtype. */
  IOStoragetype storagetype; /*!< The IOStoragetype. */

public:
  /* Sets defaults, and the valuetype matching T: */
  template <class T> IOHeader &def(const T &ref);
  IOHeader &def(const int &ref);
  IOHeader &def(const double &ref);
  IOHeader &def(const MatrixC &ref);
  IOHeader &def();
  /* Default values: */
  template <class T>
  IOHeader &dense(const Matrix<T> &M,
                  IOMatrixtype matrixt = IOMatrixtype::General);
  template <class T>
  IOHeader &sparse(const SparseMatrix<T> &M,
                   IOMatrixtype matrixt = IOMatrixtype::General);
  IOHeader &collection(const MatrixC &C);

  /* Constructor, that sets the valuetype matching T: */
  template <class T> IOHeader(const T &ref);
  IOHeader();
};

std::ostream &operator<<(std::ostream &output, const IOHeader &h);
std::istream &operator>>(std::istream &output, IOHeader &h);

/*! Base helper for input and output. */
template <class T> class IOHelper {
private:
  IOHeader h_;
  bool binary_;

public:
  /* Constructors: */
  IOHelper() : h_(T()), binary_(BINARY_DEFAULT){};
  IOHelper(const IOHeader &h) : h_(h), binary_(BINARY_DEFAULT){};

  IOHeader & header() { return h_; };
  bool is_binary() const { return binary_; };
  IOMatrixtype matrixtype() const { return h_.matrixtype; };
  IOStoragetype storagetype() const { return h_.storagetype; };

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
        (set_rowmajor ? IOStoragetype::Rowmajor : IOStoragetype::Colmajor);
    return *this;
  };
  IOHelper<T> &colmajor(bool set_colmajor = true) {
    h_.storagetype =
        (set_colmajor ? IOStoragetype::Colmajor : IOStoragetype::Rowmajor);
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
  const Matrix<T> *cM_{NULL};
  Matrix<T> *M_{NULL};

public:
  /* Constructors: */
  IOHelperM() : IOHelper<T>() {};
  IOHelperM(const IOHeader &h) : IOHelper<T>(h) {};
  IOHelperM<T> &cD(const Matrix<T> *M) {
    cM_ = M;
    M_ = NULL;
    IOHelper<T>::header().dense(*M);
    return *this;
  };
  IOHelperM<T> &D(Matrix<T> *M) {
    cM_ = M;
    M_ = M;
    IOHelper<T>::header().dense(*M);
    return *this;
  };

  IOHelperM<T> &matrixtype(IOMatrixtype matrixt) {
    IOHelper<T>::header().dense(*cM_, matrixt);
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
  IOHelperM<T> &general() { return matrixtype(IOMatrixtype::General); };
  IOHelperM<T> &symmetric() { return matrixtype(IOMatrixtype::Symmetric); };
  IOHelperM<T> &diagonal() { return matrixtype(IOMatrixtype::Diagonal); };
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
    IOHelper<T>::header().sparse(*M);
    IOHelper<T>::colmajor();
    return *this;
  };
  IOHelperSM<T> &D(SparseMatrix<T> *M) {
    cM_ = M;
    M_ = M;
    IOHelper<T>::header().sparse(*M);
    IOHelper<T>::colmajor();
    return *this;
  };

  IOHelperSM<T> &matrixtype(IOMatrixtype matrixt) {
    IOHelper<T>::header().sparse(*cM_, matrixt);
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
  IOHelperSM<T> &general() { return matrixtype(IOMatrixtype::General); };
  IOHelperSM<T> &symmetric() { return matrixtype(IOMatrixtype::Symmetric); };
  IOHelperSM<T> &diagonal() { return matrixtype(IOMatrixtype::Diagonal); };
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
    IOHelper<int>::header().collection(*M);
    return *this;
  };
  IOHelperC &D(MatrixC *M) {
    cM_ = M;
    M_ = M;
    IOHelper<int>::header().collection(*M);
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
      : loaded(false), active(false), datatype(IODatatype::Dense),
        valuetype(IOValuetype::Int), matrixtype(IOMatrixtype::General),
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
      : info(false, false, IODatatype::Dense, IOValuetype::Int,
             IOMatrixtype::General, false),
        matrix_() {
    create_blank();
  }
  MCC(IODatatype data, IOValuetype value, IOMatrixtype matrixt)
    : info(true, false, data, value, matrixt, true),
      matrix_() {
    create_blank();
  }
  template <class MatrixType>
  MCC(IODatatype data, IOValuetype value, IOMatrixtype matrixt,
      MatrixType * M,
      bool isowner = true)
    : info(true, false, data, value, matrixt, isowner), matrix_() {
    set(M, matrixt, isowner);
  }
  ~MCC() {
  }

  // Uses the existing info to create a new owned matrix
  void create_blank() {
    info.owner = true;
    if (info.datatype == IODatatype::Dense) {
      if (info.valuetype == IOValuetype::Int) {
        matrix_ = std::make_unique<Matrix<int>>();
      } else {
        matrix_ = std::make_unique<Matrix<double>>();
      }
    } else {
      if (info.valuetype == IOValuetype::Int) {
        matrix_ = std::make_unique<SparseMatrix<int>>();
      } else {
        matrix_ = std::make_unique<SparseMatrix<double>>();
      }
    }
    if (std::holds_alternative<std::unique_ptr<Matrix<int>>>(matrix_)) {
      info.datatype = IODatatype::Dense;
      info.valuetype = IOValuetype::Int;
    } else if (std::holds_alternative<std::unique_ptr<Matrix<double>>>(matrix_)) {
      info.datatype = IODatatype::Dense;
      info.valuetype = IOValuetype::Double;
    } else if (std::holds_alternative<std::unique_ptr<SparseMatrix<int>>>(matrix_)) {
      info.datatype = IODatatype::Sparse;
      info.valuetype = IOValuetype::Int;
    } else if (std::holds_alternative<std::unique_ptr<SparseMatrix<double>>>(matrix_)) {
      info.datatype = IODatatype::Sparse;
      info.valuetype = IOValuetype::Double;
    }
  }

  // Add a new matrix from raw pointer
  template <class MatrixType>
  void set(MatrixType * M,
      IOMatrixtype matrixt = IOMatrixtype::General,
      bool isowner = true) {
    info.matrixtype = matrixt;
    info.owner = isowner;
    matrix_ = M;
    if (std::holds_alternative<Matrix<int>*>(matrix_)) {
      info.datatype = IODatatype::Dense;
      info.valuetype = IOValuetype::Int;
    } else if (std::holds_alternative<Matrix<double>*>(matrix_)) {
      info.datatype = IODatatype::Dense;
      info.valuetype = IOValuetype::Double;
    } else if (std::holds_alternative<SparseMatrix<int>*>(matrix_)) {
      info.datatype = IODatatype::Sparse;
      info.valuetype = IOValuetype::Int;
    } else if (std::holds_alternative<SparseMatrix<double>*>(matrix_)) {
      info.datatype = IODatatype::Sparse;
      info.valuetype = IOValuetype::Double;
    }
    if (isowner) {
      matrix_ = std::unique_ptr<MatrixType>(M);
    }
  }

  template <class TheType>
  TheType & get() {
    if (auto ret = std::get_if<std::unique_ptr<TheType>>(&matrix_)) {
      if (ret && *ret) {
        return **ret;
      }
    } else if (auto ret = std::get_if<TheType*>(&matrix_)) {
      if (ret && *ret) {
        return **ret;
      }
    }
    set(new TheType(), info.matrixtype, true);
    return get<TheType>();
  }

  template <class TheType>
  const TheType * get_if() const {
    if (auto ret = std::get_if<std::unique_ptr<TheType>>(&matrix_)) {
      if (ret && *ret) {
        return &(**ret);
      }
    } else if (auto ret = std::get_if<TheType*>(&matrix_)) {
      if (ret && *ret) {
        return *ret;
      }
    }
    return NULL;
  }

  Matrix<int> &DI() {
    return get<Matrix<int>>();
  }
  Matrix<double> &DD() {
    return get<Matrix<double>>();
  }
  SparseMatrix<int> &SI() {
    return get<SparseMatrix<int>>();
  }
  SparseMatrix<double> &SD() {
    return get<SparseMatrix<double>>();
  }

  const Matrix<int> *cDI() const {
    return get_if<Matrix<int>>();
  }
  const Matrix<double> *cDD() const {
    return get_if<Matrix<double>>();
  }
  const SparseMatrix<int> *cSI() const {
    return get_if<SparseMatrix<int>>();
  }
  const SparseMatrix<double> *cSD() const {
    return get_if<SparseMatrix<double>>();
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
      IOMatrixtype matrixt = IOMatrixtype::General);
  template <class T>
  SparseMatrix<T> &attach(
      std::string name,
      SparseMatrix<T> *M,
      bool transfer_ownership = true,
      IOMatrixtype matrixt = IOMatrixtype::General);

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
