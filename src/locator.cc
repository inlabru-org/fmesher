/*
 *  Copyright Finn Lindgren (2010-2024)
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public License,
 *  v. 2.0. If a copy of the MPL was not distributed with this file, You can
 *  obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include <cstddef>
#include <cstring>
#include <map>
#include <set>
#include <sstream>

#include "locator.h"

namespace fmesh {

TriangleLocator::TriangleLocator(const Mesh *mesh,
                                 const std::vector<int> &dimensions,
                                 bool use_interval_tree)
  : mesh_(mesh), dim_(dimensions), bbox_(),
    bbox_locator_(dimensions.size(), use_interval_tree) {
  bbox_.resize(dim_.size());
  if (mesh_) {
    for (size_t i = 0; i < dim_.size(); ++i) {
      bbox_[i].resize(mesh_->nT());
    }

    /* Build boxes: */
    int d;
    Point mini;
    Point maxi;
    std::pair<double, double> range;
    for (size_t t = 0; t < mesh_->nT(); ++t) {
      mesh_->triangleBoundingBox(t, mini, maxi);
      for (size_t di = 0; di < dim_.size(); ++di) {
        d = dim_[di];
        range.first = mini[d];
        range.second = maxi[d];
        bbox_[di][t] = range;
      }
    }
  }

  bbox_locator_.init(bbox_.begin());
}

TriangleLocator::~TriangleLocator() { /* Nothing to do. */
}

int TriangleLocator::locate(const Point &s) const {
  FMLOG("Looking for s=" << s << std::endl);
  std::vector<double> loc(dim_.size());
  for (size_t di = 0; di < dim_.size(); ++di) {
    loc[di] = s[dim_[di]];
  }
  Dart d;
  for (bbox_locator_type::search_iterator si = bbox_locator_.search_begin(loc);
       !si.is_null(); ++si) {
    FMLOG("Starting at " << *si << std::endl);
    d = mesh_->locate_point(Dart(*mesh_, (*si)), s);
    FMLOG("Resulting dart " << d << std::endl);
    if (!d.isnull()) {
      FMLOG("Found dart " << d << std::endl);
      if (mesh_->type() == Mesh::Mtype::Sphere) {
        FMLOG("Mesh is spherical, check half-plane." << std::endl);
        const Point& s0(mesh_->S(mesh_->TV(d.t())[0]));
        const Point& s1(mesh_->S(mesh_->TV(d.t())[1]));
        const Point& s2(mesh_->S(mesh_->TV(d.t())[2]));
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
        double inner_prod = Vec::scalar(n0, s);
        double triangle_area = mesh_->triangleArea(s0, s1, s2) /
                               (mesh_->sphere_radius() * mesh_->sphere_radius());
        if (triangle_area > 2. * M_PI) {
          FMLOG("Large spherical triangle detected; point finding may be incorrect.")
        }
        if (inner_prod < -MESH_EPSILON) {
          FMLOG("Point outside spherical triangle." << std::endl);
          continue;
        }
        FMLOG("Point may be inside spherical triangle." << std::endl);
      }
      Point b;
      mesh_->barycentric(Dart(*mesh_, d.t()), s, b);
      FMLOG("Barycentric coordinates " << b << std::endl);
      if ((b[0] >= -10.0 * MESH_EPSILON) && (b[1] >= -10.0 * MESH_EPSILON) &&
          (b[2] >= -10.0 * MESH_EPSILON)) {
        FMLOG("Point found." << std::endl);
        return (d.t());
      } else {
        FMLOG("Mesh::locate_point reported incorrect finding." << std::endl);
      }
    }
  }
  FMLOG("Point not found, s=" << s << std::endl);
  return -1;
}

std::ostream &TriangleLocator::print(std::ostream &output) {
  return bbox_locator_.print(output);
}

std::ostream &operator<<(std::ostream &output, TriangleLocator &locator) {
  return locator.print(output);
}

TetraLocator::TetraLocator(const Mesh3 *mesh,
                           const std::vector<int> &dimensions,
                           bool use_interval_tree)
  : mesh_(mesh), dim_(dimensions), bbox_(),
    bbox_locator_(dimensions.size(), use_interval_tree) {
  bbox_.resize(dim_.size());
  if (mesh_) {
    for (size_t i = 0; i < dim_.size(); ++i) {
      bbox_[i].resize(mesh_->nT());
    }

    /* Build boxes: */
    int d;
    Point mini;
    Point maxi;
    std::pair<double, double> range;
    for (size_t t = 0; t < mesh_->nT(); ++t) {
      mesh_->tetraBoundingBox(t, mini, maxi);
      for (size_t di = 0; di < dim_.size(); ++di) {
        d = dim_[di];
        range.first = mini[d];
        range.second = maxi[d];
        bbox_[di][t] = range;
      }
    }
  }

  bbox_locator_.init(bbox_.begin());
}

TetraLocator::~TetraLocator() { /* Nothing to do. */
}

int TetraLocator::locate(const Point &s, Double4 &b) const {
  FMLOG("Looking for s=" << s << std::endl);
  std::vector<double> loc(dim_.size());
  for (size_t di = 0; di < dim_.size(); ++di) {
    loc[di] = s[dim_[di]];
  }
  // int search_steps = 0;
  for (bbox_locator_type::search_iterator si = bbox_locator_.search_begin(loc);
       !si.is_null(); ++si) {
    // search_steps++;
    Dart3 d(*mesh_, *si);
    FMLOG("Trying tetra " << *si << ", " << d << std::endl);
    if (!d.isnull()) {
      mesh_->barycentric(d, s, b);
      FMLOG("Barycentric coordinates " << b << std::endl);
      if ((b[0] >= -10.0 * MESH_EPSILON) && (b[1] >= -10.0 * MESH_EPSILON) &&
          (b[2] >= -10.0 * MESH_EPSILON) && (b[3] >= -10.0 * MESH_EPSILON)) {
        FMLOG("Found the containing tetra." << std::endl);
        // FMLOG("Search steps: " << search_steps << std::endl);
        return (d.t());
      }
    }
  }
  FMLOG("Point not found, s=" << s << std::endl);
  return -1;
}

std::ostream &TetraLocator::print(std::ostream &output) {
  return bbox_locator_.print(output);
}

std::ostream &operator<<(std::ostream &output, TetraLocator &locator) {
  return locator.print(output);
}



} /* namespace fmesh */
