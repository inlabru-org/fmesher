#ifndef QTOOL_H
#define QTOOL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <string>
#include <iostream>

#include "fmesher_debuglog.h"



template <typename Scalar, template<class> class Ordering=Eigen::AMDOrdering>
class QTool {
public:
  typedef Eigen::SparseMatrix<Scalar, Eigen::ColMajor> SMatrix;
  typedef Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor> MSMatrix;
  typedef Eigen::SimplicialLLT<SMatrix,
			       Eigen::Lower,
			       Ordering<typename SMatrix::StorageIndex > > SCholesky;
  typedef Eigen::TriangularView<const SMatrix, Eigen::Lower> LSMatrix;
  typedef Eigen::TriangularView<const SMatrix, Eigen::Upper> USMatrix;
protected:
  bool has_pattern_analysis_ = false;
  bool has_factorization_ = false;
  bool has_partial_inverse_ = false;
  SMatrix Q_;
  SCholesky LLt_;
  SMatrix S_;
public:
  QTool() { };

  /** Assign a new Q matrix
   *
   * \param new_pattern if \c true, any existing pattern analysis is
   * invalidated; otherwise, only the factorisation (if any) is
   * invalidated.
   */
  void Q(SMatrix& Q, bool new_pattern) {
    Q_ = Q;
    if (new_pattern) {
      invalidate_pattern_analysis();
    } else {
      invalidate_factorization();
    }
  };
  /** Assign a new Q matrix
   *
   * Any existing pattern analysis is invalidated.
   * \sa Q(SMatrix&)
   */
  void Q(const SMatrix& Q) {
    Q_ = Q;
    invalidate_pattern_analysis();
  };
  /** Assign a new Q matrix
   *
   * Any existing pattern analysis is invalidated.
   * \sa Q(MSMatrix&)
   */
  void Q(const MSMatrix& Q) {
    Q_ = Q;
    invalidate_pattern_analysis();
  };
  /** Direct Q access
   * \return A ref to the Q matrix
  */
  SMatrix& Q() {
    return Q_;
  };
  /** \return A const ref to the Q matrix
  */
  SMatrix const & Q() const {
    return Q_;
  };
  SCholesky const & LLt() {
    if (!has_factorization_) {
      factorize();
    }
    return LLt_;
  };
  /** Direct access to the partial inverse
   *
   * The partial inverse is computed if not already available.
   *
   * \return A const ref to the partial inverse of Q
  */
  SMatrix const & S() {
    if (!has_partial_inverse_) {
      partialInversion();
    }
    return S_;
  };

  /** Mark partial inverse and dependents as invalid */
  void invalidate_partial_inverse() {
    has_partial_inverse_ = false;
  };
  /** Mark factorisation and dependents as invalid */
  void invalidate_factorization() {
    has_factorization_ = false;
    invalidate_partial_inverse();
  };
  /** Mark pattern analysis and dependents as invalid */
  void invalidate_pattern_analysis() {
    has_pattern_analysis_ = false;
    invalidate_factorization();
  };

  /** Pattern analysis is prerequisite for factorization */
  void analyzePattern() {
    LLt_.analyzePattern(Q_.template selfadjointView< Eigen::Lower >());
    has_pattern_analysis_ = true;
    invalidate_factorization();
  };
  /** Factorisation requires pattern analysis, and is prerequisite for partial inverse
   *
   *  Automatically computes the pattern analysis of Q necessary
   */
  void factorize() {
    if (!has_pattern_analysis_) {
      analyzePattern();
    }
    LLt_.factorize(Q_.template selfadjointView< Eigen::Lower >());
    has_factorization_ = true;
  };

  /** Compute the log-determinant of Q
   *
   *  Automatically factorises Q if necessary
   */
  Scalar logDeterminant() {
    if (!has_factorization_) {
      factorize();
    }
    const SMatrix& L = LLt_.matrixL();
    Scalar accum = 0.0;
    for (typename SMatrix::Index i=0; i < Q_.rows(); ++i) {
      accum += std::log(L.coeff(i,i));
    }
    return 2.0 * accum;
  };

  /** Takahashi recursion
   *
   * Basic recursion:
   * For \f$i=n-1,\dots,0\f$, and \f$j=n-1,\dots,i\f$, if \f$L_{ji}\f$ is not known to be zero, let
   * \f{align*}{
   *   S_{ji} = S_{ij} &= \left(\mathbb{I}(i=j) / L_{ii} - \sum_{k=i+1}^{n-1} L_{ki} S_{kj} \right) / L_{ii}
   * \\									\
   *   &= \begin{cases}
   *         \left(1 / L_{ii} - \sum_{k=i+1}^{n-1} L_{ki} S_{kj} \right) / L_{ii}
   *         & \text{if $i=j$,} \\
   *         - \sum_{k=i+1}^{n-1} L_{ki} S_{kj} / L_{ii}
   *         & \text{otherwise.}
   *      \end{cases}
   * \f}
   *
   * Basic recursion pseudo-code:
   *
   *     for i = n-1, ..., 0
   *       for j = n-1, ..., i
   *         if (L[j,i] not known to be 0)
   *           S[j,i] = S[i,j] = (I(i==j)/L[i,i] - sum_{k=i+1}^{n-1} L[k,i] S[k,j] ) / L[i,i]
   *
   * Recursion implementation pseudo-code:
   *
   *     for i = n-1, ..., 0
   *       ScolI = reverse iterator from S[n-1,i]
   *       for j in nonzeros of column i of L, backwards, through j==i
   *         LcolI = reverse iterator from L[n-1,i]
   *         ScolJ = reverse iterator from S[n-1,j]
   *         Accumulate sum into ScolI.valueRef() of the products while k>i,
   *           decreasing LcolI and ScolJ in sync
   *         if (i == j)
   *           Add 1/L[i,i]^2 to ScolI.valueRef()
   *         else
   *           Assign ScolJ.valueRef()
   *         Decrease ScolI
   *
   */
  void partialInversion() {
    if (!has_factorization_) {
      factorize();
    }

    const SMatrix& L = LLt_.matrixL();
    SMatrix S = L.template selfadjointView< Eigen::Lower >();

    int64_t n = LLt_.rows();

    for (int64_t i=n-1; i >= 0; --i) {
      typename SMatrix::ReverseInnerIterator ScolI(S, i);
      for(typename SMatrix::ReverseInnerIterator LcolI0(L, i);
          LcolI0;
          --LcolI0) {
        //      int64_t j = LcolI0.row();
        //      FMLOG_("Processing i = " << i << ", j = " << LcolI0.row() << std::endl);
        typename SMatrix::ReverseInnerIterator LcolI(L, i);
        typename SMatrix::ReverseInnerIterator ScolJ(S, LcolI0.row());
        ScolI.valueRef() = 0.0;
        while (LcolI.row() > i) {
          while (ScolJ && (LcolI.row() < ScolJ.row())) {
            //    FMLOG_("SJ[" << ScolJ.row() << ", " << ScolJ.col() << "] move" << std::endl);
            --ScolJ;
          }
          if (ScolJ && ( LcolI.row() == ScolJ.row())) {
            //    FMLOG_("LI[" << LcolI.row() << ", " << LcolI.col() << "]"
            //              << " * "
            //              << "SJ[" << ScolJ.row() << ", " << ScolJ.col() << "] use"
            //              << std::endl);
            ScolI.valueRef() -= LcolI.value() * ScolJ.value();
            --ScolJ;
          }
          --LcolI;
        }
        if (i == LcolI0.row()) {
          //	FMLOG_("LI[" << LcolI.row() << ", " << LcolI.col() << "] use" << std::endl);
          ScolI.valueRef() += 1/LcolI.value();
          ScolI.valueRef() /= LcolI.value();
          //	FMLOG_("SI[" << ScolI.row() << ", " << ScolI.col() << "] = "
          //		  << ScolI.value()
          //		  << std::endl;
        } else {
          //	FMLOG_("LI[" << LcolI.row() << ", " << LcolI.col() << "] use" << std::endl);
          ScolI.valueRef() /= LcolI.value();
          //	FMLOG_("SI[" << ScolI.row() << ", " << ScolI.col() << "] = "
          //		  << ScolI.value()
          //		  << std::endl);
          /* This loop is for when there are nonzero S[k,j] elements
           * between the last nonzero L[k,i] element and S[i,j]
           */
          while (ScolJ.row() > i) {
            //	  FMLOG_("SJ[" << ScolJ.row() << ", " << ScolJ.col() << "] move at end"
            //		    << std::endl);
            --ScolJ;
          }
          ScolJ.valueRef() = ScolI.value();
          //	FMLOG_("SJ[" << ScolJ.row() << ", " << ScolJ.col() << "] = "
          //		  << "SI[" << ScolI.row() << ", " << ScolI.col() << "]"
          //		  << std::endl);
        }
        --ScolI;
      }
    }

    FMLOG("Partial inverse computed" << std::endl);
    FMLOG("Reordering output" << std::endl);

    // P Q P' = L L'
    // solve(P Q P') = S = solve(L L')
    // P solve(Q) P' = S
    // S_ = solve(Q) = P' S P
    S_ = S.twistedBy(LLt_.permutationP().inverse());
    FMLOG("Output reordered" << std::endl);
    has_partial_inverse_ = true;
  };







  /* L solves */
  template <class Lhs, class Rhs>
  void solveL(const Eigen::MatrixBase< Rhs > &b,
	      const Eigen::MatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::MatrixBase< Lhs >& result_ = const_cast< Eigen::MatrixBase< Lhs >& >(result);
    result_ = LLt_.matrixL().solve(b);
  }

  /* L' solves */
  template <class Lhs, class Rhs>
  void solveLt(const Eigen::MatrixBase< Rhs > &b,
	       const Eigen::MatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::MatrixBase< Lhs >& result_ = const_cast< Eigen::MatrixBase< Lhs >& >(result);
    result_ = LLt_.matrixU().solve(b);
  }


  /*
   * More useful to include the permutations:
   *
   * P Q P' = L L'
   * Q x = b
   * P Q x = P b
   * P Q P' P x = P b
   * L L' P x = P b
   * L (L' P x) = P b
   *
   * P' L x0 = b  : solvePtL
   * L' P x = x0 : solveLtP
   */

  /* P' L solves */
  template <class Lhs, class Rhs>
  void solvePtL(const Eigen::MatrixBase< Rhs > &b,
                const Eigen::MatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::MatrixBase< Lhs >& result_ = const_cast< Eigen::MatrixBase< Lhs >& >(result);
    result_ = LLt_.matrixL().solve(LLt_.permutationP() * b);
  }

  /* L' P solves */
  template <class Lhs, class Rhs>
  void solveLtP(const typename Eigen::MatrixBase< Rhs > &b,
                const typename Eigen::MatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::MatrixBase< Lhs >& result_ = const_cast< Eigen::MatrixBase< Lhs >& >(result);
    result_ = LLt_.permutationP().inverse() * LLt_.matrixU().solve(b);
  }

  /* Full P' L L' P solves */
  template <class Lhs, class Rhs>
  void solve(const Eigen::MatrixBase< Rhs > &b,
             const Eigen::MatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::MatrixBase< Lhs >& result_ = const_cast< Eigen::MatrixBase< Lhs >& >(result);
    result_.derived().resize(Q_.rows(), b.cols());
    result_ = LLt_.solve(b);
  }
  template <class Lhs, class Rhs>
  void solve(const Eigen::SparseMatrixBase< Rhs > &b,
             const Eigen::MatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::MatrixBase< Lhs >& result_ = const_cast< Eigen::MatrixBase< Lhs >& >(result);
    result_.derived().resize(Q_.rows(), b.cols());
    result_ = LLt_.solve(b);
  }
  template <class Lhs, class Rhs>
  void solve(const Eigen::SparseMatrixBase< Rhs > &b,
             const Eigen::SparseMatrixBase< Lhs > &result)
  {
    if (!has_factorization_) {
      factorize();
    }
    Eigen::SparseMatrixBase< Lhs >& result_ = const_cast< Eigen::SparseMatrixBase< Lhs >& >(result);
    result_.derived().resize(Q_.rows(), b.cols());
    result_ = LLt_.solve(b);
  }





};





#endif
