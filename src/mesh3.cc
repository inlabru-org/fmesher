/*
 *  Copyright Finn Lindgren (2010-2024)
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public License,
 *  v. 2.0. If a copy of the MPL was not distributed with this file, You can
 *  obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <map>
#include <set>
#include <sstream>
#include <array>

#include "mesh3.h"
#include "meshc.h"
#include "predicates.h"

#include "fmesher_debuglog.h"

using std::endl;

namespace fmesh {

/*
  Euler's formula for tetrahedralisations
 (Tet=Tetrahedra, T=Triangles, E=Edges, V=Vertices)
  V-E+T-Tet=1

 (100, 641, 1085, 531) -> 13
 (1000, 7360, 12971, 6474) -> 137
 (10000, 75547, 133958, 66952) -> 1459

 E \approx 12*V/2 = 6V, Tet/(2/3)*(3/4)?
 T \approx Tet*2
 Tet \approx 20*V/4=V*5, <= T*2

*/

Mesh3::Mesh3(Mtype manifold_type, size_t V_capacity, bool use_VT, bool use_TTi)
  : type_(manifold_type), use_VT_(use_VT), use_TTi_(use_TTi), TetVtx_(), TetTet_(),
    VtxTet_mapping_(), TetTeti_(), S_(), M_local_()
{
  make_M_local();
  if (V_capacity > 0) {
    TetVtx_.capacity(V_capacity * 5);
    TetTet_.capacity(V_capacity * 5);
    if (use_VT_)
      VtxTet_mapping_.reserve(V_capacity);
    if (use_TTi_)
      TetTeti_.capacity(V_capacity * 5);
    S_.capacity(V_capacity);
  }
}

Mesh3::~Mesh3() { clear(); }

Mesh3 &Mesh3::clear() {
  empty();
  return *this;
}

Mesh3 &Mesh3::empty() {
  TetVtx_.clear();
  TetTet_.clear();
  VtxTet_mapping_.clear();
  TetTeti_.clear();
  S_.clear();
  return *this;
}

Mesh3 &Mesh3::operator=(const Mesh3 &M) {
  clear();
  type_ = M.type_;
  useVT(M.use_VT_);
  useTTi(M.use_TTi_);
  S_set(M.S_);
  TV_set(M.TetVtx_);
  return *this;
}

Mesh3 &Mesh3::check_capacity(size_t nVc, size_t nTc) {
  if (nVc > S_.capacity()) {
    if (use_VT_)
      if (nVc > VtxTet_mapping_.capacity()) {
        VtxTet_mapping_.reserve(nVc);
      }
    S_.capacity(nVc);
  }

  if (nTc > TetVtx_.capacity()) {
    TetVtx_.capacity(nTc);
    TetTet_.capacity(nTc);
    if (use_TTi_)
      TetTeti_.capacity(nTc);
  }

  return *this;
}

Mesh3 &Mesh3::rebuildTT() {
  typedef Int3 Tri_Type;
  typedef std::map<Tri_Type, int> TriTet_Type;
  int t, vi, loop;
  Tri_Type T0, T1, T2;
  TriTet_Type::const_iterator Ti;
  TriTet_Type TriTet;
  TetTet_.rows(nT());
  /* Pass 1: */
  FMLOG("Pass 1" << std::endl);
  for (t = 0; t < (int)nT(); t++) {
    const Int4 &TVt = TetVtx_[t];
    Dart3 d(*this, t);
    for (loop = 0; loop < 4; ++loop) {
      FMLOG("Looking for tetra opposite to " << d << std::endl);
      const Int3 &local_TV = M_local_.TV(d.tri().t());
      T0 = Int3(TVt[ local_TV[0] ],
                TVt[ local_TV[1] ],
                   TVt[local_TV[2]]);
      for (vi = 0; vi < 3; ++vi) {
        T1 = Int3(TVt[ local_TV[(vi + 2) % 3] ],
                  TVt[ local_TV[(vi + 1) % 3] ],
                     TVt[local_TV[(vi + 0) % 3]]);
        Ti = TriTet.find(T1);
        if (Ti != TriTet.end()) {
          break;
        }
      }
      if (Ti != TriTet.end()) {
        /* Found neighbour */
        FMLOG("rebuildTT: Pass 1, Found neighbour" << std::endl);
        FMLOG("Looked for reverse of T0 = "<< T0[0] << "," << T0[1] << "," << T0[2] << std::endl);
        FMLOG("Found T1 = "<< T1[0] << "," << T1[1] << "," << T1[2] << std::endl);
        FMLOG("Found tetra = " << Ti->second << std::endl);
        FMLOG("Found mapping from = "
                << Ti->first[0] << ","
                << Ti->first[1] << ","
                << Ti->first[2] << ","
                << std::endl);
        TetTet_(t)[d.tri().t()] = Ti->second;
      } else { /* Either on boundary, or not added yet. */
        FMLOG("rebuildTT: Pass 1, On boundary, or not found yet" << std::endl);
        TetTet_(t)[d.tri().t()] = -1;
      }
      TriTet.insert(TriTet_Type::value_type(T0, t));
      FMLOG("Storing T0 = "<< T0[0] << "," << T0[1] << "," << T0[2] << std::endl);
      d.orbit3();
    }
  }

  /* Pass 2: */
  FMLOG("Pass 2" << std::endl);
  for (t = 0; t < (int)nT(); t++) {
    const Int4 &TVt = TetVtx_[t];
    Dart3 d(*this, t);
    for (loop = 0; loop < 4; ++loop) {
      FMLOG("Looking for tetra opposite to " << d << std::endl);
      if (TetTet_[t][d.tri().t()] >= 0) {
        d.orbit3();
        continue;
      }
      const Int3 &local_TV = M_local_.TV(d.tri().t());
      T0 = Int3(TVt[ local_TV[0] ],
                TVt[ local_TV[1] ],
                TVt[ local_TV[2] ]);
      const Int4 &TVt = TetVtx_[t];
      for (vi = 0; vi < 3; ++vi) {
        T1 = Int3(
          TVt[local_TV[(vi + 2) % 3]],
             TVt[local_TV[(vi + 1) % 3]],
                TVt[local_TV[(vi + 0) % 3]]
        );
        Ti = TriTet.find(T1);
        if (Ti != TriTet.end()) {
          break;
        }
      }
      if (Ti != TriTet.end()) {
        /* Found neighbour */
        FMLOG("rebuildTT: Pass 2, Found neighbour" << std::endl);
        FMLOG("Looked for reverse of T0 = "<< T0[0] << "," << T0[1] << "," << T0[2] << std::endl);
        FMLOG("Found T1 = "<< T1[0] << "," << T1[1] << "," << T1[2] << std::endl);
        FMLOG("Found tetra = " << Ti->second << std::endl);
        FMLOG("Found mapping from = "
                << Ti->first[0] << ","
                << Ti->first[1] << ","
                << Ti->first[2] << ","
                << std::endl);
        TetTet_(t)[d.tri().t()] = Ti->second;
      }
      d.orbit3();
    }
  }

  return *this;
}

Mesh3 &Mesh3::add_VT(const int v, const int t) {
  if ((use_VT_) && (v < (int)nV()) && (t < (int)nT())) {
    if (TetVtx_[t][0] == v) {
      VtxTet_mapping_[v].emplace(t, 0);
    } else if (TetVtx_[t][1] == v) {
      VtxTet_mapping_[v].emplace(t, 1);
    } else if (TetVtx_[t][2] == v) {
      VtxTet_mapping_[v].emplace(t, 2);
    } else if (TetVtx_[t][3] == v) {
      VtxTet_mapping_[v].emplace(t, 3);
    } else {
      /* Error! This should never happen! */
      FMLOG("ERROR: Vertex " << v << " not in tetra " << t << "\n");
    }
  }
  FMLOG("VT:" << VTO(v));
  check_VT_mapping_consistency();
  return *this;
}

Mesh3 &Mesh3::add_VT(const int v, const int t, const int vi) {
  if ((use_VT_) && (v < (int)nV()) && (t < (int)nT())) {
    FMLOG("Adding to VtxTet, (v, t, vi) = "
             << "(" << v << ", " << t << ", " << vi << ")" << std::endl);
    if (TetVtx_[t][vi] == v) {
      VtxTet_mapping_[v].emplace(t, vi);
  } else {
      /* Error! This should never happen! */
      FMLOG("ERROR: Vertex " << v << " doesn't match node vi=" <<
        vi << " in tetra " << t << std::endl);
    }
  }
  FMLOG("VtxTet:" << VTO(v));
  check_VT_mapping_consistency();
  return *this;
}

Mesh3 &Mesh3::remove_VT(const int v, const int t) {
  if ((use_VT_) && (v < (int)nV()) && (t < (int)nT())) {
    auto where = VtxTet_mapping_[v].find(t);
    if (where != VtxTet_mapping_[v].end())
      VtxTet_mapping_[v].erase(where);
  }
  check_VT_mapping_consistency();
  return *this;
}

Mesh3 &Mesh3::add_VT_tetra(const int t) {
  if ((use_VT_) && (t < (int)nT()) && (t >= 0)) {
    const Int4 &TVt = TetVtx_[t];
    for (int vi = 0; vi < 4; vi++) {
      add_VT(TVt[vi], t, vi);
    }
  }
  check_VT_mapping_consistency();
  return *this;
}

Mesh3 &Mesh3::remove_VT_tetra(const int t) {
  if ((use_VT_) && (t < (int)nT()) && (t >= 0)) {
    const Int4 &TVt = TetVtx_[t];
    for (int vi = 0; vi < 4; vi++) {
      remove_VT(TVt[vi], t);
    }
  }
  check_VT_mapping_consistency();
  return *this;
}

Mesh3 &Mesh3::add_VT_tetras(const int t_start) {
  if (use_VT_) {
    for (auto t = t_start; t < (int)nT(); t++) {
      add_VT_tetra(t);
      check_VT_mapping_consistency();
    }
  }
  check_VT_mapping_consistency();
  return *this;
}

Mesh3 &Mesh3::remove_VT_tetras(const int t_start) {
  if (use_VT_) {
    for (auto t = t_start; t < (int)nT(); t++) {
      remove_VT_tetra(t);
    }
  }
  check_VT_mapping_consistency();
  return *this;
}


Mesh3 &Mesh3::clear_VT(const int v) {
  if (use_VT_) {
    VtxTet_mapping_[v].clear();
  }
  return *this;
}

Mesh3 &Mesh3::reset_VT(const int v_start) {
  if (use_VT_) {
    VtxTet_mapping_.resize(nV());
    for (int v = v_start; v < (int)nV(); v++) {
      clear_VT(v);
    }
  }
  return *this;
}

Mesh3 &Mesh3::rebuild_VT() {
  if ((!use_VT_) || (!S_.capacity())) {
    VtxTet_mapping_.clear();
  } else {
    VtxTet_mapping_.clear();
    VtxTet_mapping_.reserve(S_.capacity());
    VtxTet_mapping_.resize(S_.rows());
    reset_VT(0);
    add_VT_tetras(0);
  }
  check_VT_mapping_consistency();
  return *this;
}

void Mesh3::check_VT_mapping_consistency() const {
  return;
  // if (!use_VT_)
  //   return;
  // for (int v = 0; v < (int)nV(); v++) {
  //   for (auto it = VtxTet_mapping_[v].begin(); it != VtxTet_mapping_[v].end(); it++) {
  //     if (it->first < 0 || it->first >= (int)nT()) {
  //       FMLOG_("ERROR: VtxTet_mapping_[" << v << "] contains invalid tetra "
  //                                   << it->first << std::endl);
  //     }
  //     if (it->second < 0 || it->second >= 4) {
  //       FMLOG_("ERROR: VtxTet_mapping_[" << v << "] contains invalid node index "
  //                                   << it->second << std::endl);
  //     }
  //     if (TetVtx_[it->first][it->second] != v) {
  //       FMLOG_("ERROR: VtxTet_mapping_[" << v << "] contains invalid node index "
  //                                   << it->second << " for tetra " << it->first
  //                                   << std::endl);
  //       FMLOG_(VTO(v));
  //       FMLOG_(" TetVtx[" << it->first << "] = ("
  //                << TetVtx_[it->first][0] << ", "
  //                << TetVtx_[it->first][1] << ", "
  //                << TetVtx_[it->first][2] << ", "
  //                << TetVtx_[it->first][3] << "), vi = "
  //                << it->second << std::endl);
  //     }
  //   }
  // }
  // for (int t = 0; t < (int)nT(); t++) {
  //   for (int vi = 0; vi < 4; vi++) {
  //     if (VtxTet_mapping_[TetVtx_[t][vi]].find(t) == VtxTet_mapping_[TetVtx_[t][vi]].end()) {
  //       FMLOG_("ERROR: VtxTet_mapping_[" << TetVtx_[t][vi] << "] does not contain "
  //                                   << t << std::endl);
  //       FMLOG_(VTO(TetVtx_[t][vi]));
  //       FMLOG_(" TetVtx[" << t << "] = ("
  //                     << TetVtx_[t][0] << ", "
  //                     << TetVtx_[t][1] << ", "
  //                     << TetVtx_[t][2] << ", "
  //                     << TetVtx_[t][3] << "), vi = "
  //                     << vi << std::endl);
  //     }
  //   }
  // }
}

Mesh3 &Mesh3::rebuildTTi() {
  int t, vi, t2, vi2;
  if (!use_TTi_) {
    TetTeti_.clear();
    return *this;
  }
  TetTeti_.rows(nT());
  if (!TetVtx_.capacity())
    return *this;
  TetTeti_.capacity(TetVtx_.capacity());
  for (t = 0; t < (int)nT(); t++) {
    for (vi = 0; vi < 4; vi++) {
      t2 = TetTet_[t][vi];
      if (t2 >= 0) {
        FMLOG("Connecting tetra " << t << " index " << vi << " and tetra " << t2 << std::endl);
        FMLOG(Dart3(*this, t) << std::endl);
        FMLOG(Dart3(*this, t2) << std::endl);
        for (vi2 = 0; (vi2 < 4) && (TetTet_[t2][vi2] != t); vi2++) {
        }
        FMLOG("Found index " << vi2 << std::endl);
        if (vi2 < 4) {
          TetTeti_(t)[vi] = vi2;
        } else {
          TetTeti_(t)[vi] = -1;
          /* Error! This should never happen! */
          FMLOG("ERROR\n");
        }
      } else {
        TetTeti_(t)[vi] = -1;
      }
    }
  }
  return *this;
}

Mesh3 &Mesh3::useVT(bool use_VT) {
  if (use_VT_ != use_VT) {
    use_VT_ = use_VT;
    rebuild_VT();
  }
  return *this;
}

Mesh3 &Mesh3::useTTi(bool use_TTi) {
  if (use_TTi_ != use_TTi) {
    use_TTi_ = use_TTi;
    rebuildTTi();
  }
  return *this;
}

SparseMatrix<int> Mesh3::VV() const {
  SparseMatrix<int> VV;
  for (int t = 0; t < (int)nT(); t++) {
    for (size_t i = 0; i < 4; i++) {
      for (size_t j = 0; j < 4; j++) {
        if (i != j) {
          VV(TetVtx_[t][i], TetVtx_[t][j]) = 1;
        };
      };
    };
  };
  return VV;
}

Mesh3 &Mesh3::S_set(const Matrix3double &S) {
  S_.rows(0); /* Avoid possible unnecessary copy. */
  S_append(S);
  return *this;
}

Mesh3 &Mesh3::TV_set(const Matrix4int &TV) {
  TetVtx_.rows(0); /* Avoid possible unnecessary copy. */
  TV_append(TV);
  return *this;
}

Mesh3 &Mesh3::S_append(const Point &s) {
  S_(nV()) = s;
  if (use_VT_)
    reset_VT(nV() - 1);
  return *this;
}

Mesh3 &Mesh3::S_append(const Matrix3double &S) {
  S_.append(S);
  if (use_VT_)
    reset_VT(nV() - S.rows());
  return *this;
}


Mesh3 &Mesh3::TV_append(const Matrix4int &TV) {
  TetVtx_.append(TV);
  if (use_VT_)
    add_VT_tetras(nT() - TV.rows());
  rebuildTT();
  rebuildTTi();
  return *this;
}

void Mesh3::tetraBoundingBox(const Point &s0,
                             const Point &s1,
                             const Point &s2,
                             const Point &s3,
                             Point &mini,
                             Point &maxi) const {
  for (int d = 0; d < 3; d++) {
    mini[d] = (s0[d] < s1[d] ? s0[d] : s1[d]);
    mini[d] = (mini[d] < s2[d] ? mini[d] : s2[d]);
    mini[d] = (mini[d] < s3[d] ? mini[d] : s3[d]);
    maxi[d] = (s0[d] > s1[d] ? s0[d] : s1[d]);
    maxi[d] = (maxi[d] > s2[d] ? maxi[d] : s2[d]);
    maxi[d] = (maxi[d] > s3[d] ? maxi[d] : s3[d]);
  }
}

/*!
 \brief Calculate the length of an edge.

 For planes and triangular manifolds, the edge length is
 \f$L=\|s_1-s_0\|\f$.

 On the sphere, the "obvious" arccos-formula for the length of a
 geodesic may be numerically unstable for short and long edges.  Use
 the arctan-formula instead, that should handle all cases
 \f$L\in[0,\pi R]\f$.
 \f{align*}{
 L &= R \operatorname{acos}(s_0 \cdot s_1 / R^2) \\
 \sin(L/2/R) &= \|s_1-s_0\|/2 /R \\
 \cos(L/2/R) &= \|s_0+s_1\|/2 /R \\
 L &= 2R \cdot \operatorname{atan2}(\|s_1-s_0\|,\|s_0+s_1\|)
 \f}
 */
double Mesh3::edgeLength(const Point &s0, const Point &s1) const {
  Point e;
  Vec::diff(e, s1, s0);
  double len = Vec::length(e);

  return len;
}

/*!

  \see Mesh3::edgeLength(const Dart3& d)
*/
double Mesh3::edgeLength(const Dart3 &d) const {
  int t(d.t());
  if ((t < 0) || (t >= (int)nT()))
    return 0.0;

  return edgeLength(S_[d.v()], S_[d.vo()]);
}


double Mesh3::tetraVolume(const Point &s0,
                          const Point &s1,
                          const Point &s2,
                          const Point &s3) const {
  Point e0, e1, e2;
  Vec::diff(e0, s0, s3);
  Vec::diff(e1, s1, s3);
  Vec::diff(e2, s2, s3);

  double volume = Vec::volume(e0, e1, e2) / 6.0;

  return volume;
}

void Mesh3::tetraBoundingBox(int t, Point &mini, Point &maxi) const {
  if ((t < 0) || (t >= (int)nT())) {
    return;
  }

  const auto &TVt = TetVtx_[t];
  const Point &s0 = S_[TVt[0]];
  const Point &s1 = S_[TVt[1]];
  const Point &s2 = S_[TVt[2]];
  const Point &s3 = S_[TVt[3]];
  Mesh3::tetraBoundingBox(s0, s1, s2, s3, mini, maxi);
}

double Mesh3::tetraVolume(int t) const {
  if ((t < 0) || (t >= (int)nT()))
    return 0.0;

  const auto &TVt = TetVtx_[t];
  int v0 = TVt[0];
  int v1 = TVt[1];
  int v2 = TVt[2];
  int v3 = TVt[3];
  return Mesh3::tetraVolume(S_[v0], S_[v1], S_[v2], S_[v3]);
}

void Mesh3::triangleCircumcenter(const Point &s0,
                                 const Point &s1,
                                 const Point &s2,
                                 Point &c) const {
  Point e0, e1, e2;
  Vec::diff(e0, s2, s1);
  Vec::diff(e1, s0, s2);
  Vec::diff(e2, s1, s0);

  Point n0, n1, n2;
  Vec::cross(n0, e1, e2);
  Vec::cross(n1, e2, e0);
  Vec::cross(n2, e0, e1);
  Vec::accum(n0, n1);
  Vec::accum(n0, n2);
  double scale(-4.5 / Vec::scalar(n0, n0));
  Vec::scale(c, s0, scale * Vec::scalar(e0, e0) * Vec::scalar(e1, e2));
  Vec::accum(c, s1, scale * Vec::scalar(e1, e1) * Vec::scalar(e2, e0));
  Vec::accum(c, s2, scale * Vec::scalar(e2, e2) * Vec::scalar(e0, e1));
  return;
}

/*!
 Calculate tetra circumsphere

 */
void Mesh3::tetraCircumsphere(const Point &s0,
                              const Point &s1,
                              const Point &s2,
                              const Point &s3,
                              Point &c) const {
  Point e0, e1, e2;
  Vec::diff(e0, s2, s1);
  Vec::diff(e1, s0, s2);
  Vec::diff(e2, s1, s0);

  Point n0, n1, n2;
  Vec::cross(n0, e1, e2);
  Vec::cross(n1, e2, e0);
  Vec::cross(n2, e0, e1);
  Vec::accum(n0, n1);
  Vec::accum(n0, n2);
  double scale(-4.5 / Vec::scalar(n0, n0));
  Vec::scale(c, s0, scale * Vec::scalar(e0, e0) * Vec::scalar(e1, e2));
  Vec::accum(c, s1, scale * Vec::scalar(e1, e1) * Vec::scalar(e2, e0));
  Vec::accum(c, s2, scale * Vec::scalar(e2, e2) * Vec::scalar(e0, e1));

  // We now have the circumcentre of the 012 triangle
  // Extend this in the normal direction to find the centre
  // of the circumsphere, matching the distance to s3
  // |c + a*n0 - s_3|^2 = |c + a*n0 - s_i|^2, i = 0,1,2
  //
  // |c - s_3|^2 + 2a n0.(c - s_3) + a^2 |n0|^2 =
  // |c - s_i|^2 + 2a n0.(c - s_i) + a^2 |n0|^2
  //
  // |c - s_3|^2 + 2a n0.(c - s_3) =
  // |c - s_i|^2 + 2a n0.(c - s_i)
  //
  // 2a n0.(s_avg - s_3) =
  // (|c - s_0|^2+|c - s_1|^2+|c - s_2|^2)/3 - |c - s_3|^2
  Point s_avg(s0);
  s_avg.rescale(1.0/3.0);
  s_avg.accum(s1, 1.0/3.0);
  s_avg.accum(s2, 1.0/3.0);

  Point s_tmp;
  s_tmp.diff(s_avg, s3);
  double a = Vec::scalar(n0, s_tmp);
  s_tmp.diff(c, s0);
  double b = Vec::scalar(s_tmp, s_tmp);
  s_tmp.diff(c, s1);
  b += Vec::scalar(s_tmp, s_tmp);
  s_tmp.diff(c, s2);
  b += Vec::scalar(s_tmp, s_tmp);
  s_tmp.diff(c, s3);
  b /= 3.0;
  b -= Vec::scalar(s_tmp, s_tmp);

  c.accum(n0, a / b / 2.0);
  return;
}

/*!
 Calculate tetra circumsphere

 */
void Mesh3::tetraCircumsphere(int t, Point &c) const {
  if ((t < 0) || (t >= (int)nT())) {
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 0.0;
    return;
  }

  int v0 = TetVtx_[t][0];
  int v1 = TetVtx_[t][1];
  int v2 = TetVtx_[t][2];
  int v3 = TetVtx_[t][3];
  const Point &s0 = S_[v0];
  const Point &s1 = S_[v1];
  const Point &s2 = S_[v2];
  const Point &s3 = S_[v3];

  tetraCircumsphere(s0, s1, s2, s3, c);
  return;
}

/*!
  \brief Calculate the radius of the tetra circumsphere

  We use the formula given at
  http://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates_from_cross-_and_dot-products
 \n Rewriting in our notation, the radius of the circumcircle is given by
  \f{align*}{
  r  &= \frac{3\|e_0\| \|e_1\| \|e_2\|}{2\|n_0+n_1+n_2\|}
 \f}

  \see Mesh::triangleArea
  \see Mesh::triangleCircumcenter
 */
double Mesh3::tetraCircumsphereRadius(const Point &s0,
                                      const Point &s1,
                                      const Point &s2,
                                      const Point &s3) const {
  Point c;
  Mesh3::tetraCircumsphere(s0, s1, s2, s3, c);

  Point s;
  s.diff(s3, c);
  double radius = s.length();
  return radius;
}

/*!
  \brief Calculate the radius of the tetra circumsphere
 */
double Mesh3::tetraCircumsphereRadius(int t) const {
  if ((t < 0) || (t >= (int)nT()))
    return -1.0;

  int v0 = TetVtx_[t][0];
  int v1 = TetVtx_[t][1];
  int v2 = TetVtx_[t][2];
  int v3 = TetVtx_[t][3];
  const Point &s0 = S_[v0];
  const Point &s1 = S_[v1];
  const Point &s2 = S_[v2];
  const Point &s3 = S_[v3];

  return tetraCircumsphereRadius(s0, s1, s2, s3);
}

bool Mesh3::tetraEdgeLengths(int t, Vector<double, 6> &len) const {
  if ((t < 0) || (t >= (int)nT()))
    return false;

  len[0] = edgeLength(S_[TetVtx_[t][1]], S_[TetVtx_[t][2]]);
  len[1] = edgeLength(S_[TetVtx_[t][2]], S_[TetVtx_[t][3]]);
  len[2] = edgeLength(S_[TetVtx_[t][3]], S_[TetVtx_[t][1]]);
  len[3] = edgeLength(S_[TetVtx_[t][0]], S_[TetVtx_[t][1]]);
  len[4] = edgeLength(S_[TetVtx_[t][0]], S_[TetVtx_[t][2]]);
  len[5] = edgeLength(S_[TetVtx_[t][0]], S_[TetVtx_[t][3]]);

  return true;
}

int Mesh3::tetraEdgeLengthsArgMin(int t, Vector<double, 6> &len) const {
  if (!Mesh3::tetraEdgeLengths(t, len))
    return -1;

  int amin = 0;
  for (size_t k = 0; k < 6; k++) {
    if (len[k] < len[amin]) {
      amin = k;
    }
  }
  return amin;
}

int Mesh3::tetraEdgeLengthsArgMax(int t, Vector<double, 6> &len) const {
  if (!Mesh3::tetraEdgeLengths(t, len))
    return -1;

  int amax = 0;
  for (size_t k = 0; k < 6; k++) {
    if (len[k] > len[amax]) {
      amax = k;
    }
  }
  return amax;
}


double Mesh3::tetraShortestEdge(int t) const {
  Vector<double, 6> len;
  if (!Mesh3::tetraEdgeLengths(t, len))
    return -1;

  int amin = 0;
  for (size_t k = 0; k < 6; k++) {
    if (len[k] < len[amin]) {
      amin = k;
    }
  }
  return len[amin];
}

double Mesh3::tetraLongestEdge(int t) const {
  Vector<double, 6> len;
  if (!Mesh3::tetraEdgeLengths(t, len))
    return -1;

  int amax = 0;
  for (size_t k = 0; k < 6; k++) {
    if (len[k] > len[amax]) {
      amax = k;
    }
  }
  return len[amax];
}

double Mesh3::triangleEncroached(const Dart3 &d, const Point &s) const
/* > --> encroached */
{
  int t(d.t());
  if ((t < 0) || (t >= (int)nT()))
    return -1.0;

  Dart3 dh(d);
  const Point &s0 = S_[dh.v()];
  dh.orbit2();
  const Point &s1 = S_[dh.v()];
  dh.orbit2();
  const Point &s2 = S_[dh.v()];

  Point c;
  triangleCircumcenter(s0, s1, s2, c);
  Point vec;
  vec.diff(c, s0);
  double radius = vec.length();
  Vec::diff(vec, c, s);

  /* Triangle is encroached if the distance between
   * the triangle circumcentre and the point is
   * smaller than the triangle circumcirle radius. */
  return (radius - Vec::length(vec));
}

/**!
 * Note that the tetra triangles are outwards facing.
 * Points inside the tetra are "below" the triangles.
 * Positive = inside
 */
double Mesh3::inInsideHalfspace(const Point &s0,
                                const Point &s1,
                                const Point &s2,
                                const Point &s) const {
  // TODO: check the sign of orient3d(s0, s1, s2, s3)
  return predicates::orient3d(s0.raw(),
                              s1.raw(),
                              s2.raw(),
                              s.raw());
}

double Dart3::inInsideHalfspace(const Point &s) const {
  if (isnull())
    return 0.0; /* TODO: should show a warning somewhere... */
  Dart3 dh(*this);
  int v0(dh.v());
  dh.orbit2();
  int v1(dh.v());
  dh.orbit2();
  int v2(dh.v());
  return M_->inInsideHalfspace(M_->S_[v0], M_->S_[v1], M_->S_[v2], s);
}

double Dart3::inCircumsphere(const Point &s) const {
  if (isnull())
    return 0.0; /* TODO: should show a warning somewhere... */
  const auto &TVt = M_->TetVtx_[tet_];

  return predicates::insphere(M_->S_[TVt[0]].raw(),
                              M_->S_[TVt[1]].raw(),
                              M_->S_[TVt[2]].raw(),
                              M_->S_[TVt[3]].raw(),
                              s.raw());
}



Mesh3 &Mesh3::unlinkTriangle(const Dart3 &d) {
  Dart3 dh(d);
  if (!d.onBoundary()) {
    dh.opposite3();
    TetTet_(dh.t())[dh.tri().t()] = -1;
    if (use_TTi_)
      TetTeti_(dh.t())[dh.tri().t()] = -1;
    dh = d;
  }
  TetTet_(dh.t())[dh.tri().t()] = -1;
  if (use_TTi_)
    TetTeti_(dh.t())[dh.tri().t()] = -1;

  return *this;
}

/*!
  Unlink a tetra
 */
Mesh3 &Mesh3::unlinkTetra(const int t) {
  Dart3 dh(*this, t);
  unlinkTriangle(dh);
  dh.orbit3();
  unlinkTriangle(dh);
  dh.orbit3();
  unlinkTriangle(dh);
  dh.orbit3();
  unlinkTriangle(dh);
  if (use_VT_)
    remove_VT_tetra(t);
  return *this;
}

Mesh3 &Mesh3::relocateTetra(const int t_source, const int t_target) {
  if (t_target == t_source)
    return *this;
  if (use_VT_) {
    remove_VT_tetra(t_source);
  }
  if (t_target > t_source)
    check_capacity(0, t_target + 1);
  TetVtx_(t_target)[0] = TetVtx_[t_source][0];
  TetVtx_(t_target)[1] = TetVtx_[t_source][1];
  TetVtx_(t_target)[2] = TetVtx_[t_source][2];
  TetVtx_(t_target)[3] = TetVtx_[t_source][3];
  TetTet_(t_target)[0] = TetTet_[t_source][0];
  TetTet_(t_target)[1] = TetTet_[t_source][1];
  TetTet_(t_target)[2] = TetTet_[t_source][2];
  TetTet_(t_target)[3] = TetTet_[t_source][3];
  if (use_VT_) {
    add_VT_tetra(t_target);
  }
  if (use_TTi_) {
    TetTeti_(t_target)[0] = TetTeti_[t_source][0];
    TetTeti_(t_target)[1] = TetTeti_[t_source][1];
    TetTeti_(t_target)[2] = TetTeti_[t_source][2];
    TetTeti_(t_target)[3] = TetTeti_[t_source][3];
  }
  /* Relink neighbouring TT:s. TTi is not affected by the relocation. */
  Dart3 dh(*this, t_target);
  for (size_t face = 0; face < 4; ++face) {
    if (!dh.onBoundary()) {
      dh.opposite3();
      TetTet_(dh.t())[dh.tri().t()] = t_target;
      dh.opposite3();
    }
    dh.orbit3();
  }

  return *this;
}

/*!
  Remove a tetra.

  The single-VT implementation of VT was slow when useVT is true.  For better
  performance when removing many triangles, set to false while
  removing.
  The new set-VT implementation has not been tested for speed.
 */
int Mesh3::removeTetra(const int t) {
  if ((t < 0) || (t >= (int)nT()))
    return -1;

  unlinkTetra(t);
  relocateTetra(nT() - 1, t);
  TetVtx_.rows(nT() - 1);
  /* Note: nT() was changed by the alteration of TV_. */
  TetTet_.rows(nT());
  if (use_TTi_)
    TetTeti_.rows(nT());
//  if (use_VT_)
//    rebuild_VT(); // This shouldn't be needed for the new VT implementation.
  return nT();
}

/*!
  Calculate barycentric coordinates.

*/
void Mesh3::barycentric(const Dart3 &d, const Point &s, Double4 &bary) const {
  const auto &TVt = TetVtx_[d.t()];
  int v0(TVt[0]);
  int v1(TVt[1]);
  int v2(TVt[2]);
  int v3(TVt[3]);
  bary[0] = tetraVolume(S_[v3], S_[v2], S_[v1], s);
  bary[1] = tetraVolume(S_[v2], S_[v3], S_[v0], s);
  bary[2] = tetraVolume(S_[v1], S_[v0], S_[v3], s);
  bary[3] = tetraVolume(S_[v0], S_[v1], S_[v2], s);

  Vec::rescale(bary, 1.0 / (bary[0] + bary[1] + bary[2] + bary[3]));
}

/*!
  Find the triangle opposite a vertex that a straight path will pass through.

  If the point is not found, a null Dart is returned.

  \verbatim
  findPathDirection(d0,s)
   1. d0 determines the starting vertex, v0 = d0.v
   2. d = d0
   4. d.orbit0rev
   5. while (d != d0) // Stop at edge or after one full orbit.
   6.   d.orbit0rev
   7. d.orbit2
   8. onleft0 = inLeftHalfSpace(S(v0),s,S(d.v))
   9. d.orbit2
  10. onleft1 = inLeftHalfSpace(S(v0),s,S(d.v))
  11. while !(!onleft0 & onleft1) & !d.onBoundary
  12.   d.orbit0rev
  13.   if d == d0 // We've gone a full orbit without success
  14.     return null
  15.   onleft0 = onleft1
  16.   d.orbit2
  17.   onleft1 = inLeftHalfSpace(S(v0),s,S(d.v))
  18. if (!onleft0 & onleft1)
  19.   return d
  20. return null // Not found (hit boundary).
  \endverbatim
 */
// Dart3 Mesh3::find_path_direction(const Dart &d0,
//                                  const Point &s,
//                                  const int v) const {
//   Dart3 d(d0);
//   if (d.isnull())
//     return Dart3();
//
//   int v0(d.v());
//   if (d.v() == v) // Have we found a preexisting vertex?
//     return d;
//   // TODO:
//   NOT_IMPLEMENTED;
//   d.orbit0rev();
//   while ((d != d0) && (!d.onBoundary())) {
//     d.orbit0rev();
//   }
//   Dart d00(d);
//
//   FMLOG("Finding direction to point or v starting from d00, S:"
//         << endl
//         << "\t\t" << s << endl
//         << "\t\t" << v << endl
//         << "\t\t" << d00 << endl
//         << "\t\t" << S_[d00.v()] << endl);
//
//   d.orbit2();
//   FMLOG("Finding direction to point or v starting from dart:"
//         << endl
//         << "\t\t" << s << endl
//         << "\t\t" << v << endl
//         << "\t\t" << d << endl
//         << "\t\t" << S_[d.v()] << endl
//         << "\tiLHS\t" << inLeftHalfspace(S_[v0], s, S_[d.v()]) << endl);
//   if (d.v() == v) // Have we found a preexisting vertex?
//     return d;
//   bool onleft0(inLeftHalfspace(S_[v0], s, S_[d.v()]) >= 0.0);
//   bool onleft2(d.inLeftHalfspace(s) >= -MESH_EPSILON);
//   FMLOG(d << endl);
//   d.orbit2();
//   FMLOG("Finding direction to point or v starting from dart:"
//         << endl
//         << "\t\t" << s << endl
//         << "\t\t" << v << endl
//         << "\t\t" << d << endl
//         << "\t\t" << S_[d.v()] << endl
//         << "\tiLHS\t" << inLeftHalfspace(S_[v0], s, S_[d.v()]) << endl);
//   if (d.v() == v) // Have we found a preexisting vertex?
//     return d;
//   bool onleft1(inLeftHalfspace(S_[v0], s, S_[d.v()]) >= 0.0);
//   FMLOG("Locating direction " << onleft0 << onleft1 << endl);
//   while (!(!onleft0 && onleft1) && (!d.onBoundary())) {
//     d.orbit0rev();
//     FMLOG(d << endl);
//     if (d.v() == d00.vo()) {
//       if (onleft2) {
//         FMLOG("Went full circle. Point found." << endl);
//         return d;
//       } else {
//         FMLOG("Went full circle. Point not found." << endl);
//         return Dart();
//       }
//     }
//     onleft0 = onleft1;
//     onleft2 = (onleft2 && (d.inLeftHalfspace(s) >= -MESH_EPSILON));
//     d.orbit2();
//     onleft1 = (inLeftHalfspace(S_[v0], s, S_[d.v()]) >= 0.0);
//     if (d.v() == v) // Have we found a preexisting vertex?
//       return d;
//     FMLOG("Locating direction " << onleft0 << onleft1 << endl);
//   }
//   if (!onleft0 && onleft1) {
//     d.orbit2rev();
//     return d;
//   }
//   return Dart();
// }

/*!
  Find the triangle that a straight path will pass through.

  If the path endpoint is inside the original tetra, a null Dart
  is returned.
*/
// Dart3 Mesh3::find_path_direction(const Point &s0,
//                                  const Point &s1,
//                                  const Dart3 &d0) const {
//   if (d0.isnull())
//     return Dart3();
//   Dart3 d(*this, d0.t(), 1, 0);
//
//   /* Check if we're starting on a vertex, and call alternative method */
//   /* onleft[i] = is triangle vertex i to the left of the line? */
//   /* above[i] = is s1 above triangle i? */
//   bool onleft[3];
//   for (int i = 0; i < 4; i++) {
//     if (edgeLength(S_[d.v()], s0) < MESH_EPSILON) {
//       d = find_path_direction(d, s1, -1);
//       /* Check that the line actually crosses the dart. */
//       if (inAboveHalfspace(d, s1) < 0.0) {
//         FMLOG("Checkpoint 4" << endl);
//         return d;
//       } else {
//         return Dart3();
//       }
//     }
//     // TODO
//     NOT_IMPLEMENTED;
//     onleft[i] = (inAboveHalfspace(s0, s1, S_[d.v()]) >= 0.0);
//     FMLOG("D=" << d << endl);
//     FMLOG("onleft[" << i << "] = " << onleft[i] << endl);
//     d.orbit2();
//   }
//
//   // TODO
//   NOT_IMPLEMENTED;
//   for (int i0 = 0; i0 < 3; i0++) {
//     int i1 = (i0 + 1) % 3;
//     if (inLeftHalfspace(S_[d.v()], S_[d.vo()], s1) < 0.0) {
//       if (!onleft[i1]) {
//         FMLOG("Checkpoint 1" << endl);
//         d.orbit2();
//         return d;
//       }
//       if (onleft[i0]) {
//         FMLOG("Checkpoint 2" << endl);
//         d.orbit2rev();
//         return d;
//       }
//       FMLOG("Checkpoint 3" << endl);
//       return d;
//     }
//     d.orbit2();
//   }
//   return Dart3();
// }

/*!
  Trace the geodesic path from a vertex to a point or another vertex.

  Return a pair of darts identifying the starting vertex, together
  with the point containing triangle, or a dart originating at the
  endpoint vertex.  If the point/vertex is not found, a null Dart is
  returned.  Priority is given to finding the vertex; if found, the
  point is disregarded.  The trace only includes darts strictly
  intersected, i.e. not the initial and final darts.

  Alg 9.1 is non-robust, and does not take the shortest route.  This
  algorithm finds and follows the straight-line intersected edges
  instead, making it suitable for use in CDT segment insertion
  algorithms.

  \verbatim
  tracePath:
   1. d0 determines the starting vertex, v0 = d0.v
   2. d = findPathDirection(d0,s) // d is now opposite v0
   3. if inLeftHalfspace(d,s) return d
   4. while !d.onBoundary
   5.   store d in path-trace
   6.   d1 = d.orbit1.orbit2rev
   7.   found1 = inLeftHalfspace(d1,s)
   8.   leavethrough2 = inLeftHalfspace(S(v0),s,S(d.v))
   9.   d2 = d1.orbit2rev
  10.   found2 = inLeftHalfspace(d2,s)
  11.   if found1 & found2 return d2
  12.   if leavethrough2, d=d2, else d=d1
  13. store d in path-trace
  14. return null
  \endverbatim
 */
// Dart3Pair Mesh3::trace_path(const Dart3 &d0,
//                             const Point &s1,
//                             const int v1,
//                             Dart3List *trace) const {
//   Dart3 dh;
//   if (d0.isnull())
//     dh = Dart3(*this, 0);
//   else
//     dh = Dart3(*this, d0.t(), 1, d0.vi());
//   int v0(dh.v());
//   FMLOG("Locating point " << s1 << " v0=" << v0 << " v1=" << v1 << endl);
//
//   if (v1 >= (int)nV()) { /* Vertex index out of range */
//     return Dart3Pair(dh, Dart3());
//   }
//
//   Dart3 d(find_path_direction(dh, s1, v1));
//   FMLOG("Path-direction " << d << endl);
//   FMLOG("Starting triangle "
//           << d.t() << " ("
//           << TetVtx_[d.t()][0] << ","
//           << TetVtx_[d.t()][1] << ","
//           << TetVtx_[d.t()][2] << ")"
//           << TetVtx_[d.t()][3] << ")"
//           << endl);
//   if (d.isnull()) {
//     FMLOG("No direction found, so is in starting triangle" << endl);
//     return Dart3Pair(dh, dh);
//   }
//   Dart3 dstart = d;
//   while (dstart.v() != d0.v())
//     dstart.orbit2rev();
//   FMLOG("Starting dart " << dstart << endl);
//   if ((d.v() == v1) || (d.inAboveHalfspace(s1) >= -MESH_EPSILON)) {
//     FMLOG("Found " << d << endl);
//     return Dart3Pair(dstart, d);
//   }
//   while (!d.onBoundary()) {
//     if (trace) {
//       trace->push_back(d);
//     }
//     // TODO:
//     NOT_IMPLEMENTED;
//     d.orbit1().orbit2rev();
//     FMLOG("In triangle " << d << endl);
//     if (d.v() == v1) {
//       FMLOG("Found vertex at " << d << endl);
//       return Dart3Pair(dstart, d);
//     }
//     bool found = (d.inAboveHalfspace(s1) >= -MESH_EPSILON);
//     bool other = (inAboveHalfspace(S_[v0], s1, S_[d.v()]) > 0.0);
//     d.orbit2rev();
//     if (found && (d.inAboveHalfspace(s1) >= -MESH_EPSILON)) {
//       return Dart3Pair(dstart, d);
//     }
//     if (!other)
//       d.orbit2();
//     FMLOG("Go to next triangle, from " << d << endl);
//   }
//   FMLOG("Endpoint not found " << dstart << " " << d << endl);
//   return Dart3Pair(dstart, Dart3());
// }

/*!
  Trace the geodesic path between two points.

  Return a pair of darts identifying the starting and end triangles.
  If the end point is not reached due to reaching a mesh boundary
  the second dart returned will be null, but the path trace up to an
  including the mesh boundary will be populated.  The trace only
  includes darts strictly intersected.

  \verbatim
  tracePath:
   2. d = findPathDirection(s0,s1,d0)
   3. if d is null (i.e. inLeftHalfspace(d,s1)), return (d0,d0)
   4. while !d.onBoundary
   5.   store d in path-trace
   6.   d1 = d.orbit1.orbit2rev
   7.   found1 = inLeftHalfspace(d1,s1)
   8.   leavethrough2 = inLeftHalfspace(s0,s1,S(d.v))
   9.   d2 = d1.orbit2rev
  10.   found2 = inLeftHalfspace(d2,s1)
  11.   if found1 & found2 return (d0,d2)
  12.   if leavethrough2, d=d2, else d=d1
  13. store d in path-trace
  14. return (d0,null)
  \endverbatim
 */
// Dart3Pair Mesh3::trace_path(const Point &s0, const Point &s1, const Dart3 &d0,
//                           Dart3List *trace) const {
//   // TODO:
//   NOT_IMPLEMENTED;
//   Dart dh;
//   if (d0.isnull()) {
//     return DartPair(Dart(), Dart());
//   }
//   dh = Dart(*this, d0.t(), 1, d0.vi());
//   FMLOG("Locating point s1=" << s1 << "from s0=" << s0 << endl);
//   Dart dstart = dh;
//   FMLOG("Starting dart " << dstart << endl);
//
//   Dart d(find_path_direction(s0, s1, dstart));
//   FMLOG("Path-direction " << d << endl);
//   FMLOG("Starting triangle " << d.t() << " (" << TV_[d.t()][0] << ","
//                              << TV_[d.t()][1] << "," << TV_[d.t()][2] << ")"
//                              << endl);
//   if (d.isnull()) {
//     FMLOG("No direction found, so is in starting triangle" << endl);
//     return DartPair(dstart, dstart);
//   }
//   if (d.inLeftHalfspace(s1) >= -MESH_EPSILON) {
//     FMLOG("Found " << d << endl);
//     return DartPair(dstart, dstart);
//   }
//   while (!d.onBoundary()) {
//     if (trace) {
//       trace->push_back(d);
//     }
//     d.orbit1().orbit2rev();
//     FMLOG("In triangle " << d << endl);
//     bool found = (d.inLeftHalfspace(s1) >= -MESH_EPSILON);
//     bool other = (inLeftHalfspace(s0, s1, S_[d.v()]) > 0.0);
//     d.orbit2rev();
//     if (found && (d.inLeftHalfspace(s1) >= -MESH_EPSILON)) {
//       return DartPair(dstart, d);
//     }
//     if (!other)
//       d.orbit2();
//     FMLOG("Go to next triangle, from " << d << endl);
//   }
//   FMLOG("Endpoint not found " << dstart << " " << d << endl);
//   if (trace) {
//     trace->push_back(d);
//   }
//   return DartPair(dstart, Dart());
// }

/*!
  Locate a point in the graph.

  Return a dart identifying the containing tetra.

  If the point is not found, a null Dart3 is returned.
 */
// Dart3 Mesh3::locate_point(const Dart3 &d0, const Point &s, const int v) const {
//   Dart3 dh;
//   if (d0.isnull()) {
//     dh = Dart3(*this, 0); /* Another option here would be to pick a
//                             random tetra instead of always tetra 0.
//                             For now, we leave such options to the caller.
//                          */
//   } else {
//     dh = Dart3(*this, d0.t(), 1, d0.vi());
//   }
//   return trace_path(dh, s, v).second;
// }

/*!
  Locate an existing vertex in the graph.

  Return a dart originating at the vertex.

  If the vertex is not found, a null Dart3 is returned.
 */
Dart3 Mesh3::locate_vertex(const Dart3 &d0, const int v) const {
  if ((v < 0) || (v >= (int)nV())) {
    return Dart3(); /* Vertex index out of range */
  }

  if (use_VT_) {
    if (VtxTet_mapping_[v].empty()) {
      /* Vertex not connected to any triangles. */
      return Dart3();
    }

    /* Return arbitrary dart from the vertex. */
    auto v_t_iter = VtxTet_mapping_[v].begin();
    Dart d_local(M_local_, v_t_iter->second, 1, 0);
    // d_local is now in the triangle _not_ containing the vertex!
    // Move to the neighbouring triangle, and then the opposite vertex.
    d_local.orbit1().orbit2rev();
    if (d_local.v() != v_t_iter->second) {
      FMLOG_("Error: Vertex not found in tetra containing it." << std::endl);
    }
    return Dart3(*this, v_t_iter->first, d_local);
  }

  NOT_IMPLEMENTED;
  FMLOG_("use_VT_ must be true in Mesh3::locate_vertex()"
           << std::endl);
  Dart3 d(d0);
  if (d.isnull())
    return Dart3();
  return Dart3();

  // Dart3 dh;
  // if (d0.isnull())
  //   dh = Dart3(*this, 0);
  // else
  //   dh = Dart3(*this, d0.t(), 1, d0.vi());
  // dh = trace_path(dh, S_[v], v).second;
  // if (dh.v() != v) /* Point may be found, but not the actual vertex. */
  //   return Dart3();
  // return dh;
}

MOAint4 Mesh3::TVO() const { return MOAint4(TetVtx_, nT()); }
MOAint4 Mesh3::TTO() const { return MOAint4(TetTet_, nT()); }
MOAVTMap Mesh3::VTO() const { return MOAVTMap(VtxTet_mapping_, nV()); }
MOAVTMapV Mesh3::VTO(const int v) const { return MOAVTMapV(VtxTet_mapping_, v); }
MOAint4 Mesh3::TTiO() const { return MOAint4(TetTeti_, nT()); }
MOAdouble3 Mesh3::SO() const { return MOAdouble3(S_, nV()); }

void Mesh3::calcQblocks(SparseMatrix<double> &C0, SparseMatrix<double> &C1,
                        SparseMatrix<double> &G1, SparseMatrix<double> &B1,
                        Matrix<double> &Tvolumes) const {
  C0.clear().rows(nV()).cols(nV());
  C1.clear().rows(nV()).cols(nV());
  G1.clear().rows(nV()).cols(nV());
  B1.clear().rows(nV()).cols(nV());
  Tvolumes.clear().cols(1).rows(nT());
  for (int t = 0; t < (int)nT(); t++) {
    const Int4Raw &tv = TetVtx_[t].raw();
    const Point &s0 = S_[tv[0]];
    const Point &s1 = S_[tv[1]];
    const Point &s2 = S_[tv[2]];
    const Point &s3 = S_[tv[3]];
    Point e[4];
    e[0].diff(s3, s2);
    e[1].diff(s1, s3);
    e[2].diff(s2, s1);
    e[3].diff(s0, s1);

    // TODO: should likely be pairwise product of triangle normals
    NOT_IMPLEMENTED;
    std::array<std::array<double, 4>, 4> eij;
    for (size_t i = 0; i < 4; i++) {
      eij[i][i] = Vec::scalar(e[i], e[i]);
      for (size_t j = i + 1; j < 4; j++) {
        eij[i][j] = Vec::scalar(e[i], e[j]);
        eij[j][i] = eij[i][j];
      }
    }

    bool b[4];
    b[0] = (TetTet_[t][0] < 0 ? true : false);
    b[1] = (TetTet_[t][1] < 0 ? true : false);
    b[2] = (TetTet_[t][2] < 0 ? true : false);
    b[3] = (TetTet_[t][3] < 0 ? true : false);

    double a = tetraVolume(t);
    Tvolumes(t, 0) = a;

    // Flat volume; no curved 3D manifolds supported.
    double fa = a;

    double vij;
    for (int i = 0; i < 4; i++) {
      C0(tv[i], tv[i]) += a / 4.;
      C1(tv[i], tv[i]) += a / 6.; // TODO
      G1(tv[i], tv[i]) += eij[i][i] / (4. * fa); // TODO
      for (int j = i + 1; j < 4; j++) {
        C1(tv[i], tv[j]) += a / 12.; // TODO
        C1(tv[j], tv[i]) += a / 12.; // TODO
        vij = eij[i][j] / (4. * fa); // TODO
        G1(tv[i], tv[j]) += vij; // TODO
        G1(tv[j], tv[i]) += vij; // TODO
      }
    }

    // TODO
    NOT_IMPLEMENTED;
    if (b[0] || b[1] || b[2]) {
      vij = -1. / (4. * fa);
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            if (b[k] && (i != k)) {
              B1(tv[i], tv[j]) += eij[k][j] * vij;
            }
          }
        }
      }
    }
  }
}

void crossmultiply3(const Vector<double, 4> *ax, Vector<double, 4> *H, int n) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      H[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        H[i][j] += ax[k][i] * ax[k][j];
      }
    }
  }
}

void adjugate3(const Vector<double, 4> *H, Vector<double, 4> *aH) {
  // TODO
  NOT_IMPLEMENTED;
  aH[0][0] = H[1][1] * H[2][2] - H[1][2] * H[2][1];
  aH[0][1] = H[1][2] * H[2][0] - H[1][0] * H[2][2];
  aH[0][2] = H[1][0] * H[2][1] - H[1][1] * H[2][0];
  aH[1][1] = H[0][0] * H[2][2] - H[0][2] * H[2][0];
  aH[1][2] = H[0][1] * H[2][0] - H[0][0] * H[2][1];
  aH[2][2] = H[0][0] * H[1][1] - H[0][1] * H[1][0];
  aH[1][0] = aH[0][1];
  aH[2][0] = aH[0][2];
  aH[2][1] = aH[1][2];
}



// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R

bool Mesh3::save(std::string filename_s, std::string filename_tv,
                bool binary) const {
  return (S_.save(filename_s, IOMatrixtype::General, binary) &&
          TetVtx_.save(filename_tv, IOMatrixtype::General, binary));
}

bool Mesh3::load(std::string filename_s, std::string filename_tv, bool binary) {
  Matrix3<double> s;
  Matrix4<int> tv;

  if (!s.load(filename_s, binary) || !tv.load(filename_tv, binary))
    return false;

  S_set(s);
  TV_set(tv);

  return true;
}

#endif // not FMESHER_WITH_R

// Change the vertex only
// Same tetra
// Same triangle
// Same edge but opposite direction
Dart3 &Dart3::alpha0() {
  tri_.alpha0();
  return *this;
}

// Same vertex
// Change the edge only (will reverse orientation)
// Same triangle
// Same tetra
Dart3 &Dart3::alpha1() {
  tri_.alpha1();
  return *this;
}

// Same vertex
// Same edge (will reverse orientation)
// Change the triangle only
// Same tetra
Dart3 &Dart3::alpha2() {
  tri_.alpha2();
  return *this;
}

Dart3 &Dart3::alpha3() {
  if (onBoundary())
    return *this;
  opposite3().alpha1();
  return *this;
}

// Same vertex
// Same triangle (but seen from other side)
// Different edge
// Opposite tetra
Dart3 &Dart3::opposite3() {
  if (onBoundary())
    return *this;
  int tet0 = tet_;
  int v0 = M_->TetVtx_[tet_][tri_.v()];
  int tet = M_->TetTet_[tet_][tri_.t()];
  if (tet < 0)
    return *this;
  int vi;
  // Seek the opposing vertex/triangle
  if (!M_->use_TTi_) {
    for (vi = 0; (vi < 4) && (M_->TetTet_[tet][vi] != tet0); vi++) {
    }
    if (vi >= 4) {
      return *this; /* Error! This should never happen! */
    }
  } else {
    vi = M_->TetTeti_[tet0][tri_.t()];
  }
  tet_ = tet;
  tri_ = Dart(M_->M_local_, vi, tri_.edir(), 0);
  // Seek the appropriate triangle rotation
  if (M_->TetVtx_[tet][tri_.v()] != v0) {
    tri_.orbit2();
    if (M_->TetVtx_[tet][tri_.v()] != v0) {
      tri_.orbit2();
      if (M_->TetVtx_[tet][tri_.v()] != v0) {
        FMLOG_("Error! This should never happen!" << std::endl);
        return *this; /* Error! This should never happen! */
      }
    }
  }

  return *this;
}

// Go to the next triangle in the same tetra
// with the vertex in common
Dart3 &Dart3::orbit0() {
  tri_.orbit0();
  return *this;
}

// Go to the adjacent triangle in the same tetra
// with a common edge.
Dart3 &Dart3::orbit1() {
  tri_.orbit1();
  return *this;
}

// Go to the next point in the triangle.
Dart3 &Dart3::orbit2() {
  tri_.orbit2();
  return *this;
}

// Move around in the tetra.
Dart3 &Dart3::orbit3() {
  tri_ = Dart(M_->M_local_, (tri_.t() + 1) % 4, tri_.edir(), tri_.vi());
  return *this;
}

Dart3 &Dart3::orbit0rev() {
  tri_.orbit0rev();
  return *this;
}

Dart3 &Dart3::orbit1rev() /* Equivalent to orbit1() */
{
  tri_.orbit1rev();
  return *this;
}

Dart3 &Dart3::orbit2rev() {
  tri_.orbit2rev();
  return *this;
}

// Move around in the tetra.
Dart3 &Dart3::orbit3rev() {
  tri_ = Dart(M_->M_local_, (tri_.t() + 4 - 1) % 4, tri_.edir(), tri_.vi());
  return *this;
}

std::ostream &operator<<(std::ostream &output, const Mesh3 &M) {
  output << "Mesh type:\t" << M.type() << endl;
  output << "Vertices:\t" << M.nV() << endl;
  output << "Tetrahedra:\t" << M.nT() << endl;
  output << "Options:\t" << (M.useVT() ? "VT " : "")
         << (M.useTTi() ? "TTi " : "")
  << endl;
  return output;
}

std::ostream &operator<<(std::ostream &output, const Mesh3::Mtype &type) {
  switch (type) {
  case Mesh3::Mtype::Manifold:
    output << "Manifold (R3)";
    break;
  case Mesh3::Mtype::Plane:
    output << "Plane (R3)";
    break;
  }
  return output;
}

std::ostream &operator<<(std::ostream &output, const MOAint4 &MO) {
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < (int)MO.n_; i++) {
      output << ' ' << std::right << std::setw(4) << MO.M_[i][j];
    }
    output << endl;
  }
  return output;
}

std::ostream &operator<<(std::ostream &output, const MOAdouble4 &MO) {
  for (int i = 0; i < (int)MO.n_; i++) {
    for (int j = 0; j < 4; j++)
      output << ' ' << std::right << std::setw(10) << std::scientific
             << MO.M_[i][j];
    output << endl;
  }
  return output;
}

std::ostream &operator<<(std::ostream &output, const Vector<double, 4> &MO) {
  output << '(';
  for (int j = 0; j < 4; j++) {
    output << std::right << std::setw(10) << std::scientific << MO[j];
    if (j < 2)
      output << ',';
  }
  output << ')';
  return output;
}

std::ostream &operator<<(std::ostream &output, const Dart3 &d) {
  output << "D3=(" << std::right << std::setw(1) << d.tet_
         << ", " << d.tri_
         << ")";
  if ((!d.isnull()) && (d.tet_ < (int)d.M()->nT())) {
    output << " TetV=("
           << d.M()->TV(d.tet_)[0] << ","
           << d.M()->TV(d.tet_)[1] << ","
           << d.M()->TV(d.tet_)[2] << ","
           << d.M()->TV(d.tet_)[3] << ")";
    output << " TT=("
           << d.M()->TT(d.tet_)[0] << ","
           << d.M()->TT(d.tet_)[1] << ","
           << d.M()->TT(d.tet_)[2] << ","
           << d.M()->TT(d.tet_)[3] << ")";
  }

  return output;
}

} /* namespace fmesh */
