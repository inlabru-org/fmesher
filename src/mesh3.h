/*
 *  Copyright Finn Lindgren (2010-2024)
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public License,
 *  v. 2.0. If a copy of the MPL was not distributed with this file, You can
 *  obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef _FMESH_MESH3_
#define _FMESH_MESH3_ 1

#include <cstddef>
//#include <cstring>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "fmesher_debuglog.h"
#include "vector.h"
#include "mesh.h"

#define MESH3_EPSILON 1e-15

namespace fmesh {

class Mesh3;
class Dart3;
class MOAint4;
class MOAdouble4;

// typedef std::pair<int, int> IntPair;
// typedef std::list<int> vertexListT;
typedef std::set<int> tetraSetT;
typedef std::list<Dart3> Dart3List;
typedef std::pair<Dart3, Dart3> Dart3Pair;
// typedef std::map<int,int> IntMap;
// typedef std::vector<IntMap> VTMapT;

class Mesh3 {
  friend class Dart3;
  friend std::ostream &operator<<(std::ostream &output, const Mesh3 &M);

public:
  enum class Mtype : int { Manifold = 0, Plane};
  friend std::ostream &operator<<(std::ostream &output, const Mesh3::Mtype &type);

private:
  Mtype type_;
  bool use_VT_;
  bool use_TTi_;
  Matrix4int TetVtx_;  /* TetVtx[t]  : {v0,v1,v2,v3} */
  Matrix4int TetTet_;  /* TT[t]  : {t0,t1,t2,t3} */
  VTMapT VtxTet_mapping_;  /* VtxTet[v] : map from t to the vi (0,1,2,3) value for v */
  /* TetVtx[ VtxTet[v][i].t, VtxTet[v][i].vi ] == v */
  Matrix4int TetTeti_; /* TetTeti[t] : {vi0,vi1,vi2,vi3},
                     t == TetTet[ TetTet[t][i] ][ TetTeti[t][i] ] */
  Matrix3double S_;
  Mesh M_local_;

private:
  Mesh3 &rebuildTT();

  /*! Add tetrahedron to VtxTet[v] */
  Mesh3 &add_VT(const int v, const int t);
  /*! Add tetra to VtxTet[v] with precomputed vi information */
  Mesh3 &add_VT(const int v, const int t, const int vi);
  /*! Remove tetra to VtxTet[v] */
  Mesh3 &remove_VT(const int v, const int t);
  /*! Clear VtxTet[v] info */
  Mesh3 &clear_VT(const int v);
  /*! Set VtxTet[v]=empty for v>=v_start */
  Mesh3 &reset_VT(const int v_start = 0);

  Mesh3 &add_VT_tetra(const int t);
  Mesh3 &remove_VT_tetra(const int t);
  Mesh3 &add_VT_tetras(const int t_start = 0);
  Mesh3 &remove_VT_tetras(const int t_start = 0);
  void check_VT_mapping_consistency() const;

  Mesh3 &rebuild_VT();
  Mesh3 &rebuildTTi();

  void make_M_local() {
    Matrix3double S_local;
    Matrix3int TV_local;
    S_local.rows(4);
    TV_local.rows(4);
    S_local(0) = Point(1.0, 0.0, 0.0);
    S_local(1) = Point(0.0, 1.0, 0.0);
    S_local(2) = Point(0.0, 0.0, 1.0);
    S_local(3) = Point(0.0, 0.0, 0.0);
    TV_local(0) = Int3(3, 2, 1);
    TV_local(1) = Int3(2, 3, 0);
    TV_local(2) = Int3(1, 0, 3);
    TV_local(3) = Int3(0, 1, 2);

    M_local_.type(Mesh::Mtype::Manifold);
    M_local_.S_set(S_local);
    M_local_.TV_set(TV_local);
    M_local_.useVT(true);
    M_local_.useTTi(true);
  }

public:
  Mesh3(void)
    : type_(Mtype::Manifold), use_VT_(false), use_TTi_(true),
      TetVtx_(), TetTet_(),
      VtxTet_mapping_(), TetTeti_(), S_(), M_local_()
  {
    make_M_local();
  };
  Mesh3(Mtype manifold_type, size_t Vcapacity, bool use_VT = true,
       bool use_TTi = false);
  Mesh3(const Mesh3 &M)
      : type_(Mtype::Manifold), use_VT_(true),
        use_TTi_(false), TetVtx_(), TetTet_(), VtxTet_mapping_(), TetTeti_(), S_(),
        M_local_()
  {
    *this = M;
  };
  Mesh3 &operator=(const Mesh3 &M);
  ~Mesh3();
  Mesh3 &clear();
  Mesh3 &empty();

  /*!
    \brief Check the storage capacity, and increase if necessary
  */
  Mesh3 &check_capacity(size_t nVc, size_t nTc);
  size_t Vcap() const { return S_.capacity(); }

  bool useVT() const { return use_VT_; };
  Mesh3 &useVT(bool use_VT);
  bool useTTi() const { return use_TTi_; };
  Mesh3 &useTTi(bool use_TTi);

  Mtype type() const { return type_; };
  void type(Mtype set_type) { type_ = set_type; };

  Mtype determine_type() {
    return Mtype::Plane;
  }


  Mtype auto_type() {
    Mtype mtype = determine_type();
    (*this).type(mtype);
    return mtype;
  }




  size_t nV() const { return S_.rows(); };
  size_t nT() const { return TetVtx_.rows(); };
  const Matrix4int &TV() const { return TetVtx_; };
  const Matrix4int &TT() const { return TetTet_; };
  const VTMapT &VT() const { return VtxTet_mapping_; };
  const Matrix4int &TTi() const { return TetTeti_; };
  const Matrix3double &S() const { return S_; };

  Matrix4int &TV() { return TetVtx_; };
  Matrix4int &TT() { return TetTet_; };
  VTMapT &VT() { return VtxTet_mapping_; };
  Matrix4int &TTi() { return TetTeti_; };
  Matrix3double &S() { return S_; };

  SparseMatrix<int> VV() const;
  const Int4 &TV(int t) const { return TetVtx_[t]; };
  const Int4 &TT(int t) const { return TetTet_[t]; };
  const IntMap &VT(int v) const { return VtxTet_mapping_[v]; };
  const Int4 &TTi(int t) const { return TetTeti_[t]; };
  const Point &S(int v) const { return S_[v]; };

  MOAint4 TVO() const;
  MOAint4 TTO() const;
  MOAVTMap VTO() const;
  MOAVTMapV VTO(const int v) const;
  MOAint4 TTiO() const;
  MOAdouble3 SO() const;

  Mesh3 &S_set(const Matrix3double &S);
  Mesh3 &TV_set(const Matrix4int &TV);
  Mesh3 &S_append(const Point &s);
  Mesh3 &S_append(const Matrix3double &S);
  Mesh3 &TV_append(const Matrix4int &TV);

  /*
   Dart3 find_path_direction(const Dart3 &d0,
                            const Point &s,
                            const int v = -1) const;
  Dart3 find_path_direction(const Point &s0,
                            const Point &s1,
                            const Dart3 &d0) const;
  Dart3Pair trace_path(const Dart3 &d0,
                       const Point &s,
                       const int v = -1,
                       Dart3List *trace = NULL) const;
  Dart3Pair trace_path(const Point &s0,
                       const Point &s1,
                       const Dart3 &d0,
                       Dart3List *trace = NULL) const;
  Dart3 locate_point(const Dart3 &d0,
                     const Point &s,
                     const int v = -1) const;
   */
  Dart3 locate_vertex(const Dart3 &d0, const int v) const;

  Mesh3 &removeLastVertex() { /* Does not check that the vertex is unused! */
    if (nV() > 0) {
      S_.rows(nV() - 1);
      if (use_VT_)
        VtxTet_mapping_.pop_back();
    }
    return *this;
  };
  Mesh3 &unlinkTriangle(const Dart3 &d);
  Mesh3 &unlinkTetra(const int t);
  Mesh3 &relocateTetra(const int t_source,
                       const int t_target);
  int removeTetra(const int t);

  /* Traits: */
  double edgeLength(const Point &s0,
                    const Point &s1) const;
  double edgeLength(const Dart3 &d) const;
  void barycentric(const Dart3 &d,
                   const Point &s,
                   Double4 &bary) const;
  void tetraBoundingBox(const Point &s0,
                        const Point &s1,
                        const Point &s2,
                        const Point &s3,
                        Point &mini, Point &maxi) const;
  double tetraVolume(const Point &s0,
                     const Point &s1,
                     const Point &s2,
                     const Point &s3) const;
  void tetraCircumsphere(const Point &s0,
                         const Point &s1,
                         const Point &s2,
                         const Point &s3,
                         Point &c) const;
  void triangleCircumcenter(const Point &s0,
                            const Point &s1,
                            const Point &s2,
                            Point &c) const;
  double tetraCircumsphereRadius(const Point &s0,
                                 const Point &s1,
                                 const Point &s2,
                                 const Point &s3) const;
  double edgeIntersection(const Point &s00,
                          const Point &s01,
                          const Point &s10,
                          const Point &s11,
                          Point &c) const;

  void tetraBoundingBox(int t, Point &mini, Point &maxi) const;
  double tetraVolume(int t) const;
  void tetraCircumsphere(int t, Point &c) const;
  double tetraCircumsphereRadius(int t) const;
  bool tetraEdgeLengths(int t, Vector<double, 6> &len) const;
  int tetraEdgeLengthsArgMin(int t, Vector<double, 6> &len) const;
  int tetraEdgeLengthsArgMax(int t, Vector<double, 6> &len) const;
  double tetraLongestEdge(int t) const;
  double tetraShortestEdge(int t) const;
  double triangleEncroached(const Dart3 &d, const Point &s) const;

  /*!
    \brief Compute dart half-space test for a point.

    Positive if s is inside the triangle defined by d.
   */
  double inInsideHalfspace(const Point &s0,
                           const Point &s1,
                           const Point &s2,
                           const Point &s) const;

  /*!
    \brief Calculate FEM matrices.
   */
  void calcQblocks(SparseMatrix<double> &C0,
                   SparseMatrix<double> &C1,
                   SparseMatrix<double> &G1,
                   SparseMatrix<double> &B1,
                   Matrix<double> &Tvolumes) const;
  std::vector<SparseMatrix<double>> calcGradientMatrices() const;

  // No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
  /*! \brief Store the mesh in files. */
  bool save(std::string filename_s, std::string filename_tv,
            bool binary = true) const;
  /*! \brief Read a mesh from files. */
  bool load(std::string filename_s, std::string filename_tv,
            bool binary = true);
#endif // not FMESHER_WITH_R
};

class MOAint4 {
  friend std::ostream &operator<<(std::ostream &output, const MOAint4 &MO);

private:
  size_t n_;
  const Matrix4int &M_;

public:
  MOAint4(const Matrix4int &M, size_t n) : n_(n), M_(M){};
};

class MOAdouble4 {
  friend std::ostream &operator<<(std::ostream &output, const MOAdouble4 &MO);

private:
  size_t n_;
  const Matrix4double &M_;

public:
  MOAdouble4(const Matrix4double &M, size_t n) : n_(n), M_(M){};
};

/*! \brief Darts */
class Dart3 {
  friend std::ostream &operator<<(std::ostream &output, const Dart3 &d);

private:
  const Mesh3 *M_;
  Dart tri_; // Dart for a triangle, using local indexing in M_->M_local_
  // tri_.t_ is the local index of the opposing vertex in the tetra
  int tet_;

public:
  Dart3(void) : M_(NULL), tri_(), tet_(0){};
  Dart3(const Mesh3 &M, int t = 0)
    : M_(&M), tri_(), tet_(t){
    tri_ = Dart(M_->M_local_);
  };
  Dart3(const Mesh3 &M, int t, Dart tri)
    : M_(&M), tri_(tri), tet_(t){};
  Dart3(const Dart3 &d) : M_(d.M_), tri_(d.tri_), tet_(d.tet_){};
  Dart3 &operator=(const Dart3 &d) {
    M_ = d.M_;
    tri_ = d.tri_;
    tet_ = d.tet_;
    return *this;
  };

  const Mesh3 *M() const { return M_; };
  Dart tri() const { return tri_; };
  int t() const { return tet_; };
  int v() const {
    if (!M_)
      return -1;
    else
      return M_->TetVtx_[tet_][tri_.v()];
  };
  /* Vertex at other end of edge; alpha0().v() */
  int vo() const {
    if (!M_)
      return -1;
    else
      return M_->TetVtx_[tet_][tri_.vo()];
  };
  /* Adjacent tetra; opposite3().t() */
  int tadj() const {
    if (!M_)
      return -1;
    else
      return M_->TetTet_[tet_][tri_.t()];
  };

  bool isnull() const { return (!M_) || tri_.isnull(); };
  bool operator==(const Dart3 &d) const {
    return ((d.tet_ == tet_) && (d.tri_ == tri_));
  };
  bool operator<(const Dart3 &d) const {
    /* TODO: Add debug check for M_==d.M_ */
    return ((d.tet_ < tet_) ||
            ((d.tet_ == tet_) &&
             (d.tri_ < tri_)));
  };
  bool operator!=(const Dart3 &d) const { return !(d == *this); };
  bool operator>(const Dart3 &d) const { return (d < *this); };

  bool onBoundary() const { return (M_->TetTet_[tet_][tri_.t()] < 0); }

  double inInsideHalfspace(const Point &s) const;
  double inCircumsphere(const Point &s) const;
  bool circumsphereOK(void) const;

  /* Graph traversal algebra. */
  Dart3 &alpha0(void); //! Change vertex only
  Dart3 &alpha1(void); //! Change edge only
  Dart3 &alpha2(void); //! Change triangle only
  Dart3 &alpha3(void); //! Change tetra only
  Dart3 &orbit0(void); //! Orbit around vertex, fixed tetra
  Dart3 &orbit1(void); //! Orbit around edge, fixed tetra
  Dart3 &orbit2(void); //! Orbit around triangle, fixed tetra
  Dart3 &orbit3(void); //! Orbit around tetra, fixed tetra
  Dart3 &opposite3(void); //! Go to neighbouring tetra, same triangle and vertex
  Dart3 &orbit0rev(void);
  Dart3 &orbit1rev(void);
  Dart3 &orbit2rev(void);
  Dart3 &orbit3rev(void);
};

} /* namespace fmesh */

#endif
