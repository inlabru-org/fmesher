
#include <cmath>
#include <cstddef>
#include <cstring>
#include <map>
#include <set>
#include <sstream>

#include "fmesher_debuglog.h"
#include "ioutils.h"
#include "vector.h"

#include "RcppFmesher.h"

using std::endl;
using std::ios;

namespace fmesh {

// /*! dense/sparse/map */
// enum class IODatatype : int {
//   Invalid = -1,
//     Dense = 0,
//     Sparse = 1,
//     Collection = 2
// };
// /*! int/double */
// enum class IOValuetype : int {
//   Invalid = -1,
//     Int = 0,
//     Double = 1
// };
// /*! rowmajor/colmajor */
// enum class IOStoragetype : int {
//   Invalid = -1,
//     Rowmajor = 0,
//     Colmajor = 1
// };


// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R


IOHeader::IOHeader() { def(); }

IOHeader &IOHeader::def(const int &ref) {
  (void)(ref);
  def();
  valuetype = IOValuetype::Int;
  return *this;
}

IOHeader &IOHeader::def(const double &ref) {
  (void)(ref);
  def();
  valuetype = IOValuetype::Double;
  return *this;
}

IOHeader &IOHeader::def() {
  version = IOHEADER_VERSION;
  elems = 0;
  rows = 0;
  cols = 0;
  datatype = IODatatype::Invalid;
  valuetype = IOValuetype::Invalid;
  matrixtype = IOMatrixtype::Invalid;
  storagetype = IOStoragetype::Rowmajor;
  return *this;
}

IOHeader &IOHeader::collection(const MatrixC &C) {
  datatype = IODatatype::Collection;
  elems = C.output_size();
  rows = -1;
  cols = -1;
  storagetype = IOStoragetype::Invalid;
  return *this;
}

std::ostream &operator<<(std::ostream &output, const IOHeader &h) {
  int header_length = sizeof(IOHeader) / sizeof(int);
  output << header_length << " ";
  const int *ioheader_p = (const int *)&h;
  for (int i = 0; i < header_length; i++) {
    output << ioheader_p[i];
    if (i + 1 < header_length)
      output << " ";
  }
  return output;
}

std::istream &operator>>(std::istream &input, IOHeader &h) {
  int header_length = sizeof(IOHeader) / sizeof(int);
  int file_header_length;
  input >> file_header_length;
  if (file_header_length < header_length) {
    h.dense(Matrix<int>(0), IOMatrixtype::General);
    int *ioheader_p = (int *)&h;
    for (int i = 0; i < file_header_length; i++)
      input >> ioheader_p[i];
  } else {
    int *ioheader_p = (int *)&h;
    for (int i = 0; i < header_length; i++)
      input >> ioheader_p[i];
    if (file_header_length > header_length) {
      int *buf = new int[header_length - file_header_length];
      for (int i = 0; i < header_length - file_header_length; i++)
        input >> buf[i];
      delete[] buf;
    }
  }
  return input;
}

IOHelperC &IOHelperC::OL(std::ostream &output) {
  const IOHeader &h(IOHelper<int>::header());
  if ((h.elems == 0) || (!cM_)) {
    return *this;
  }
  for (const auto& outi : cM_->output_) {
    if (is_binary()) {
      int string_size = outi.length() + 1;
      output.write((char *)&string_size, sizeof(string_size));
      output.write((char *)outi.c_str(), sizeof(char) * string_size);
    } else {
      output << outi << std::endl;
    }
  }
  return *this;
}

IOHelperC &IOHelperC::IL(std::istream &input) {
  const IOHeader &h(IOHelper<int>::header());
  const bool bin_(IOHelper<int>::is_binary());
  if ((h.elems == 0) || (!cM_)) {
    return *this;
  }
  std::string name;
  for (int i = 0; i < h.elems; i++) {
    if (bin_) {
      int string_size;
      input.read((char *)&string_size, sizeof(string_size));
      char *buf = new char[string_size];
      input.read(buf, sizeof(char) * string_size);
      list_.push_back(std::string(buf));
      delete[] buf;
    } else {
      input >> name;
      list_.push_back(name);
    }
  }
  return *this;
}

IOHelperC &IOHelperC::OD(std::ostream &output) {
  const IOHeader &h(IOHelper<int>::header());
  const bool bin_(IOHelper<int>::is_binary());
  if ((h.elems == 0) || (!cM_)) {
    return *this;
  }
  for (const auto& outi : cM_->output_) {
    const MCC &mcc = *(cM_->coll_.find(outi)->second);
    if (mcc.info.datatype == IODatatype::Dense)
      if (mcc.info.valuetype == IOValuetype::Int) {
        IOHelperM<int> ioh;
        ioh.cD(mcc.cDI()).matrixtype(mcc.info.matrixtype);
        ioh.binary(bin_).OH(output).OD(output);
      } else {
        IOHelperM<double> ioh;
        ioh.cD(mcc.cDD()).matrixtype(mcc.info.matrixtype);
        ioh.binary(bin_).OH(output).OD(output);
      }
    else if (mcc.info.valuetype == IOValuetype::Int) {
      IOHelperSM<int> ioh;
      ioh.cD(mcc.cSI()).matrixtype(mcc.info.matrixtype);
      ioh.binary(bin_).OH(output).OD(output);
    } else {
      IOHelperSM<double> ioh;
      ioh.cD(mcc.cSD()).matrixtype(mcc.info.matrixtype);
      ioh.binary(bin_).OH(output).OD(output);
    }
  }
  return *this;
}

IOHelperC &IOHelperC::ID(std::istream &input) {
  const IOHeader &h(IOHelper<int>::header());
  const bool bin_(IOHelper<int>::is_binary());
  if ((!M_) || (h.elems == 0)) {
    return *this;
  }
  for (const auto& listi : list_) {

    IOHelper<int> ioh_;
    ioh_.binary(bin_).IH(input);

    if (ioh_.header().datatype == IODatatype::Dense)
      if (ioh_.header().valuetype == IOValuetype::Int) {
        IOHelperM<int> ioh;
        ioh.D(&(M_->DI(listi)));
        ioh.binary(bin_).IH(ioh_.header()).ID(input);
      } else {
        IOHelperM<double> ioh;
        ioh.D(&(M_->DD(listi)));
        ioh.binary(bin_).IH(ioh_.header()).ID(input);
      }
    else if (ioh_.header().valuetype == IOValuetype::Int) {
      IOHelperSM<int> ioh;
      ioh.D(&(M_->SI(listi)));
      ioh.binary(bin_).IH(ioh_.header()).ID(input);
    } else {
      IOHelperSM<double> ioh;
      ioh.D(&(M_->SD(listi)));
      ioh.binary(bin_).IH(ioh_.header()).ID(input);
    }
  }
  return *this;
}

#endif // not FMESHER_WITH_R



template <>
Matrix<int> &MatrixC::attach(std::string name, Matrix<int> *M,
                             IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Dense, IOValuetype::Int,
                                               matrixt, M, false)));
  activate(name);
  return coll_[name]->DI();
}


template <>
Matrix<double> &MatrixC::attach(std::string name, Matrix<double> *M,
                                IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Dense, IOValuetype::Double,
                                               matrixt, M, false)));
  activate(name);
  return coll_[name]->DD();
}


template <>
SparseMatrix<int> &MatrixC::attach(std::string name, SparseMatrix<int> *M,
                                   IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Sparse, IOValuetype::Int,
                                               matrixt, M, false)));
  activate(name);
  return coll_[name]->SI();
}

template <>
SparseMatrix<double> &MatrixC::attach(std::string name, SparseMatrix<double> *M,
                                      IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Sparse, IOValuetype::Double,
                                               matrixt, M, false)));
  activate(name);
  return coll_[name]->SD();
}

template <>
Matrix<int> &MatrixC::attach(std::string name,
                             std::unique_ptr<Matrix<int>>&& M,
                             IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Dense, IOValuetype::Int,
                                               matrixt, std::move(M))));
  activate(name);
  return coll_[name]->DI();
}

template <>
Matrix<double> &MatrixC::attach(std::string name,
                                std::unique_ptr<Matrix<double>>&& M,
                                IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Dense, IOValuetype::Double,
                                               matrixt, std::move(M))));
  activate(name);
  return coll_[name]->DD();
}

template <>
SparseMatrix<int> &MatrixC::attach(std::string name,
                                   std::unique_ptr<SparseMatrix<int>>&& M,
                                   IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Sparse, IOValuetype::Int,
                                               matrixt, std::move(M))));
  activate(name);
  return coll_[name]->SI();
}

template <>
SparseMatrix<double> &MatrixC::attach(std::string name,
                                      std::unique_ptr<SparseMatrix<double>>&& M,
                                      IOMatrixtype matrixt) {
  free(name);
  coll_.insert(collPairT(name,
                         std::make_unique<MCC>(IODatatype::Sparse, IOValuetype::Double,
                                               matrixt, std::move(M))));
  activate(name);
  return coll_[name]->SD();
}

bool MatrixC::activate(std::string name) {
  collT::iterator colli;
  if ((colli = coll_.find(name)) == coll_.end()) {
    return false;
  }
  colli->second->info.active = true;
  return true;
}

void MatrixC::activate() {
  for (const auto& colli : coll_) {
    colli.second->info.active = true;
  }
}

// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R

void MatrixC::load_file(std::string filename, bool only_list) {
  IOHelperC ioh;
  if (filename == "-") {
    /* Can only read stdin once, so read everything now. */
    ioh.D(this);
    ioh.binary(bin_in_).IH(FM_CIN).IL(FM_CIN).ID(FM_CIN);
  } else {
    std::ifstream I;
    I.open(filename.c_str(), (bin_in_ ? (ios::in | ios::binary) : ios::in));
    if (!I.is_open()) {
      // TODO: Add error handling.
    }
    ioh.D(this);
    ioh.binary(bin_in_).IH(I).IL(I);
    if (!only_list)
      ioh.ID(I);
    I.close();
  }

  /* Populate source_ */
  for (const auto& listi : ioh.list_) {
    source_[listi] = filename;
  }
}

MCCInfo MatrixC::load(std::string name) {
  /* Is the matrix already loaded? */
  if (activate(name))
    return info(name);

  sourceT::const_iterator sourcei;
  if ((sourcei = source_.find(name)) != source_.end()) {
    /* The matrix is in a collection file */
    load_file(sourcei->second);
    if (activate(name))
      return info(name);
  }

  /* Do we have a prefix to read from? */
  if (input_prefix_ == "-")
    return info(name);

  /* Try to read from prefix data. */

  std::ifstream I;
  I.open((input_prefix_ + name).c_str(),
         (bin_in_ ? (ios::in | ios::binary) : ios::in));
  if (!I.is_open()) {
    return info(name);
  }
  IOHelper<int> ioh_;
  ioh_.binary(bin_in_).IH(I);
  I.close();

  coll_.insert(collPairT(name, new MCC(IODatatype(ioh_.header().datatype),
                                       IOValuetype(ioh_.header().valuetype),
                                       IOMatrixtype(ioh_.header().matrixtype))));
  activate(name);

  I.open((input_prefix_ + name).c_str(),
         (bin_in_ ? (ios::in | ios::binary) : ios::in));
  if (!I.is_open()) {
    return info(name);
  }
  if (ioh_.header().datatype == IODatatype::Dense)
    if (ioh_.header().valuetype == IOValuetype::Int) {
      IOHelperM<int> ioh;
      ioh.D(&DI(name));
      ioh.binary(bin_in_).IH(I);
      ioh.ID(I);
    } else {
      IOHelperM<double> ioh;
      ioh.D(&DD(name));
      ioh.binary(bin_in_).IH(I).ID(I);
    }
  else if (ioh_.header().valuetype == IOValuetype::Int) {
    IOHelperSM<int> ioh;
    ioh.D(&SI(name));
    ioh.binary(bin_in_).IH(I).ID(I);
  } else {
    IOHelperSM<double> ioh;
    ioh.D(&SD(name));
    ioh.binary(bin_in_).IH(I).ID(I);
  }
  I.close();

  return info(name);
}

#endif // not FMESHER_WITH_R

MatrixC &MatrixC::free(std::string name) {
  dont_output(name);

  collT::iterator colli;
  if ((colli = coll_.find(name)) != coll_.end()) {
    coll_.erase(colli);
  }
  return *this;
}

MatrixC &MatrixC::dont_output(std::string name) {
  outputT::iterator outi;
  if ((outi = output_.find(name)) != output_.end()) {
    output_.erase(outi);
  }
  return *this;
}

MatrixC &MatrixC::output(std::string name) {
  if (name == "-") {
    output_all_ = true;
    for (auto const & colli : coll_) {
      if (colli.second->info.active)
        output_.insert(colli.first);
    }
  } else if (name == "--") {
    output_all_ = true;
    for (auto const & colli : coll_) {
      if (activate(colli.first))
        output_.insert(colli.first);
    }
  } else {
    if (info(name).loaded) {
      activate(name);
      if (output_all_) {
        output_all_ = false;
        output_.clear();
      }
      output_.insert(name);
    }
  }
  return *this;
}

void MatrixC::io(bool bin_in, bool bin_out) {
  bin_in_ = bin_in;
  bin_out_ = bin_out;
}

void MatrixC::input_prefix(std::string prefix) { input_prefix_ = prefix; }

void MatrixC::output_prefix(std::string prefix) { output_prefix_ = prefix; }

// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
void MatrixC::input_file(std::string filename) { load_file(filename, true); }
void MatrixC::output_file(std::string filename) { output_file_ = filename; }
#endif

// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
void MatrixC::input_raw(std::string name, std::string specification,
                        std::string filename) {
  /* Parse raw ascii matrix data and add to collection. */

  if (specification == "ddgr") {
    Matrix<double> &M = attach(name, new Matrix<double>());
    if (filename == "-")
      input_raw_M(std::cin, M);
    else {
      std::ifstream I;
      I.open(filename.c_str(), std::ios::in);
      if (!I.is_open()) {
        // TODO: Add error handling.
      }
      input_raw_M(I, M);
      I.close();
    }
  } else if (specification == "digr") {
    Matrix<int> &M = attach(name, new Matrix<int>());
    if (filename == "-")
      input_raw_M(std::cin, M);
    else {
      std::ifstream I;
      I.open(filename.c_str(), std::ios::in);
      if (!I.is_open()) {
        // TODO: Add error handling.
      }
      input_raw_M(I, M);
      I.close();
    }
  } else
    NOT_IMPLEMENTED;
}

void MatrixC::save() {
  /* Write the matrix collection to output */
  if (output_prefix_ != "-") {
    for (auto const & outi : output_) {
      MCC &mcc = *(coll_.find(outi)->second);
      if (mcc.info.datatype == IODatatype::Dense)
        if (mcc.info.valuetype == IOValuetype::Int)
          save_M((output_prefix_ + outi), mcc.DI(), mcc.info, bin_out_);
        else
          save_M((output_prefix_ + outi), mcc.DD(), mcc.info, bin_out_);
      else if (mcc.info.valuetype == IOValuetype::Int)
        save_SM((output_prefix_ + outi), mcc.SI(), mcc.info, bin_out_);
      else
        save_SM((output_prefix_ + outi), mcc.SD(), mcc.info, bin_out_);
    }
  }
  if (output_file_ != "") {
    if (output_file_ == "-") {
      IOHelperC ioh;
      ioh.cD(this);
      ioh.binary(bin_out_).OH(FM_COUT).OL(FM_COUT).OD(FM_COUT);
    } else {
      std::ofstream O;
      O.open(output_file_.c_str(),
             (bin_out_ ? (ios::out | ios::binary) : ios::out));
      if (!O.is_open()) {
        // TODO: Add error handling.
      }
      IOHelperC ioh;
      ioh.cD(this);
      ioh.binary(bin_out_).OH(O).OL(O).OD(O);
      O.close();
    }
  }
}
#endif // not FMESHER_WITH_R

Matrix<int> &MatrixC::DI(std::string name) {
  collT::iterator colli;
  if (((colli = coll_.find(name)) != coll_.end()) &&
      (colli->second->info.datatype == IODatatype::Dense) &&
      (colli->second->info.valuetype == IOValuetype::Int) &&
      (colli->second->info.active)) {
    return colli->second->DI();
  }
  return attach(name, std::make_unique<Matrix<int>>());
}

Matrix<double> &MatrixC::DD(std::string name) {
  collT::iterator colli;
  if (((colli = coll_.find(name)) != coll_.end()) &&
      (colli->second->info.datatype == IODatatype::Dense) &&
      (colli->second->info.valuetype == IOValuetype::Double) &&
      (colli->second->info.active)) {
    return colli->second->DD();
  }
  return attach(name, std::make_unique<Matrix<double>>());
}

SparseMatrix<int> &MatrixC::SI(std::string name) {
  collT::iterator colli;
  if (((colli = coll_.find(name)) != coll_.end()) &&
      (colli->second->info.datatype == IODatatype::Sparse) &&
      (colli->second->info.valuetype == IOValuetype::Int) &&
      (colli->second->info.active)) {
    return colli->second->SI();
  }
  return attach(name, std::make_unique<SparseMatrix<int>>());
}

SparseMatrix<double> &MatrixC::SD(std::string name) {
  collT::iterator colli;
  if (((colli = coll_.find(name)) != coll_.end()) &&
      (colli->second->info.datatype == IODatatype::Sparse) &&
      (colli->second->info.valuetype == IOValuetype::Double) &&
      (colli->second->info.active)) {
    return colli->second->SD();
  }
  return attach(name, std::make_unique<SparseMatrix<double>>());
}

void MatrixC::matrixtype(std::string name, IOMatrixtype matrixt) {
  collT::iterator colli;
  if ((colli = coll_.find(name)) != coll_.end())
    colli->second->info.matrixtype = matrixt;
}

MCCInfo MatrixC::info(std::string name) const {
  collT::const_iterator colli;
  if ((colli = coll_.find(name)) == coll_.end()) {
    return MCCInfo();
  }
  return colli->second->info;
}


#ifdef FMESHER_WITH_R
SEXP MatrixC::Rcpp_wrap() const {
  /* Convert the matrix collection to a list of R objects */
  Rcpp::List res;
  for (auto const & outi : output_) {
    const MCC &mcc = *(coll_.find(outi)->second);
    if (mcc.info.datatype == IODatatype::Dense) {
      if (mcc.info.valuetype == IOValuetype::Int) {
        if (auto ptr = mcc.cDI()) {
          res[outi] = Rcpp::wrap(*ptr);
        }
      } else {
        if (auto ptr = mcc.cDD()) {
          res[outi] = Rcpp::wrap(*ptr);
        }
      }
    } else if (mcc.info.valuetype == IOValuetype::Int) {
      if (auto ptr = mcc.cSI()) {
        res[outi] = Rcpp::wrap(*ptr);
      }
    } else {
      if (auto ptr = mcc.cSD()) {
        res[outi] = Rcpp::wrap(*ptr);
      }
    }
  }
  return res;
}


void MatrixC::attach(std::string name, SEXP from) {
  if (Rcpp::is<Rcpp::NumericMatrix>(from)) {
    (*this).attach(name, std::make_unique<Matrix<double>>(
        Rcpp::as<Rcpp::NumericMatrix>(
          from)),
          IOMatrixtype::General);
  } else if (Rcpp::is<Rcpp::IntegerMatrix>(from)) {
    (*this).attach(name, std::make_unique<Matrix<int>>(
        Rcpp::as<Rcpp::IntegerMatrix>(
          from)),
          IOMatrixtype::General);
  } else if (Rcpp::is<Rcpp::CharacterMatrix>(from)) {
  } else if (Rcpp::is<Rcpp::NumericVector>(from)) {
    (*this).attach(name, std::make_unique<Matrix<double>>(
        Rcpp::as<Rcpp::NumericVector>(
          from)),
          IOMatrixtype::General);
  } else if (Rcpp::is<Rcpp::IntegerVector>(from)) {
    (*this).attach(name, std::make_unique<Matrix<int>>(
        Rcpp::as<Rcpp::IntegerVector>(
          from)),
          IOMatrixtype::General);
  } else if (Rcpp::is<Rcpp::CharacterVector>(from)) {
  } else {
    (*this).attach(name,
     std::make_unique<SparseMatrix<double>>(from),
     IOMatrixtype::General);
  }
}




void MatrixC::attach(SEXP from) {
    Rcpp::List from_list = Rcpp::as<Rcpp::List>(from);
  Rcpp::CharacterVector from_names = from_list.names();
  for (auto const & elem : from_names) {
    std::string the_name = Rcpp::as<std::string>(elem);
    (*this).attach(the_name, from_list[the_name]);
  }
}

MatrixC::MatrixC(SEXP from)
  : output_all_(false), bin_in_(true), bin_out_(true), input_prefix_("-"),
    output_prefix_("-"), output_file_("") {
  (*this).attach(from);
}



#endif

} /* namespace fmesh */
