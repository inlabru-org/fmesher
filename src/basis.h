#ifndef _FMESH_BASIS_
#define _FMESH_BASIS_ 1

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "vector.h"

namespace fmesh {

std::unique_ptr<Matrix<double>> spherical_harmonics(
    const Matrix3<double> &S, size_t max_order,
    bool rotationally_symmetric);

std::unique_ptr<Matrix<double>> spherical_bsplines1(
    const Matrix<double> &S, size_t n_basis,
    size_t degree,
    bool uniform_knot_angle_spacing);
std::unique_ptr<Matrix<double>> spherical_bsplines(
    const Matrix3<double> &S, size_t n_basis,
    size_t degree,
    bool uniform_knot_angle_spacing);

} /* namespace fmesh */

#endif
