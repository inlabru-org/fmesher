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
#include <vector>
#ifdef FMESHER_WITH_GSL
  #include "gsl/gsl_sf_legendre.h"
#endif

#include "fmesher_debuglog.h"
#include "ioutils.h"
#include "vector.h"

#define M_2_SQRT_PI 3.5449077018110320546

using std::endl;

namespace fmesh {

int sph_basis_n(int kmax, bool rot_inv) {
  if (kmax >= 0) {
    if (rot_inv)
      return (kmax + 1);
    else
      return (kmax + 1) * (kmax + 1);
  } else {
    return 0;
  }
}

int sph_basis_index(int order, int mode, bool rot_inv) {
  return order * (order + 1) + mode;
}

size_t legendre_array_index(int order, int mode, bool rot_inv) {
#ifdef FMESHER_WITH_GSL
  return gsl_sf_legendre_array_index(order, mode);
#else
  if (rot_inv) {
    return order;
  }
  /* Ordering: order by order and mode within order */
  /* order: 0, 1,1, 2,2,2, 3,3,3,3, ... */
  /* mode:  0, 0,1, 0,1,2, 0,1,2,3, ... */
  /* index: 0, 1,2, 3,4,5, 6,7,8,9, ... */
  /* 0, 1, 3, 6, 10, ... = o*(o+1)/2 */
  return (order * (order + 1)) / 2 + mode;
#endif
}
size_t legendre_array_n(int max_order, bool rot_inv) {
#ifdef FMESHER_WITH_GSL
  return gsl_sf_legendre_array_n(max_order);
#else
  if (rot_inv) {
    return max_order + 1;
  }
  return ((max_order + 2) * (max_order + 1)) / 2;
#endif
}

/* The GSL SPHARM normalisation includes \sqrt{1/(4\pi)}
 * Using the same normalisation, in our own implementation, even though we
 * then undo the 4\pi factor in spherical_harmonics()
 */
void legendre_array(int max_order, double x,
                    double* res_array) {
#ifdef FMESHER_WITH_GSL
  gsl_sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, max_order, x, res_array);
#else
  /* P_0^0(x) = 1 / \sqrt(4\pi) */
  if (max_order <= 0) {
    /* Scaling for spherical harmonics, like GSL_SF_LEGENDRE_SPHARM */
    res_array[legendre_array_index(0, 0, true)] = 1.0 / M_2_SQRT_PI;
    return;
  }
  res_array[legendre_array_index(0, 0, true)] = 1.0;
  /* P_1^0(x) = x */
  res_array[legendre_array_index(1, 0, true)] = x;
  for (int n = 2; n <= max_order; n++) {
    /* P_n^0(x) = ((2n-1)x P_{n-1}^0(x) - (n-1) P_{n-2}^0(x)) / n */
    res_array[legendre_array_index(n, 0, true)] =
        ((2.0 * n - 1.0) * x * res_array[legendre_array_index(n - 1, 0, true)] -
         (n - 1.0) * res_array[legendre_array_index(n - 2, 0, true)]) /
        n;
  }

  /* Scaling for spherical harmonics, like GSL_SF_LEGENDRE_SPHARM */
  for (int n = 0; n <= max_order; n++) {
    res_array[legendre_array_index(n, 0, true)] *=
      std::sqrt(2.0 * n + 1.0) / M_2_SQRT_PI;
  }

#endif
}
/* Condon-Shortley phase multiplier (-1)^m when csphase = -1 */
void legendre_array_e(int max_order, double x, int csphase,
                      double* res_array) {
#ifdef FMESHER_WITH_GSL
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, max_order, x, csphase, res_array);
#else
  int max_n(legendre_array_n(max_order, false));
  for (int n = 0; n < max_n; n++) {
    res_array[n] = 0.0;
  }

  /* P_0^0(x) = 1 / \sqrt(4\pi) */
  if (max_order <= 0) {
    /* Scaling for spherical harmonics, like GSL_SF_LEGENDRE_SPHARM */
    res_array[legendre_array_index(0, 0, false)] = 1.0 / M_2_SQRT_PI;
    return;
  }
  res_array[legendre_array_index(0, 0, false)] = 1.0;
  for (int m = 0; m <= max_order; m++) {
    /* P_m^m(x) = -(2m-1) sqrt(1-x^2) P_{m-1}^{m-1}(x) */
    if (m >= 1) {
      res_array[legendre_array_index(m, m, false)] =
          csphase * (2.0 * m - 1.0) * std::sqrt(1.0 - x * x) *
          res_array[legendre_array_index(m - 1, m - 1, false)];
    }
    if (m < max_order) {
      /* P_{m+1}^m(x) = (2m+1)x P_m^m(x) */
      res_array[legendre_array_index(m + 1, m, false)] =
        (2.0 * m + 1.0) * x *
        res_array[legendre_array_index(m, m, false)];
      for (int n = m + 2; n <= max_order; n++) {
        /* P_n^m(x) = ((2n-1)x P_{n-1}^m(x) - (n+m-1) P_{n-2}^m(x)) / (n - m) */
        res_array[legendre_array_index(n, m, false)] =
        ((2.0 * n - 1.0) * x *
        res_array[legendre_array_index(n - 1, m, false)] -
        (n + m - 1.0) * res_array[legendre_array_index(n - 2, m, false)]) /
          (n - m);
      }
    }
  }

  /* Scaling for spherical harmonics, like GSL_SF_LEGENDRE_SPHARM */
  for (int n = 0; n <= max_order; n++) {
    double scaling = (2.0 * n + 1.0) / M_2_SQRT_PI / M_2_SQRT_PI;
    for (int m = 0; m <= n; m++) {
      if (m > 0) {
        scaling *= 1 / double(n + m) / double(n - m + 1);
      }
      res_array[legendre_array_index(n, m, false)] *= std::sqrt(scaling);
    }
  }

#endif
}


std::unique_ptr<Matrix<double>> spherical_harmonics(
    const Matrix3<double> &S,
    size_t max_order,
    bool rot_inv) {

  auto sph =
      std::make_unique<Matrix<double>>(sph_basis_n(max_order, rot_inv));

  size_t i, k, m;
  size_t GSL_res_n = legendre_array_n(max_order, rot_inv);
  auto GSL_res_array = std::make_unique<double[]>(GSL_res_n);

  if (rot_inv) {
    for (i = 0; i < S.rows(); i++) {
      legendre_array(max_order, S[i][2], &GSL_res_array[0]);
      for (k = 0; k <= max_order; k++) {
        (*sph)(i, k) =
            M_2_SQRT_PI * GSL_res_array[legendre_array_index(k, 0, rot_inv)];
      }
    }
  } else {
    double phi, scaling_sin, scaling_cos;

    std::vector<size_t> Idxs2(max_order + 1);
    // 0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132, 156, 182, 210, 240
    for (k = 0; k <= max_order; k++) {
      Idxs2[k] = k * (k + 1);
    }

    for (i = 0; i < S.rows(); i++) {
      phi = atan2(S[i][1], S[i][0]);

      legendre_array_e(max_order, S[i][2], -1, &GSL_res_array[0]);
      for (k = 0; k <= max_order; k++) {
        (*sph)(i, Idxs2[k]) =
            M_2_SQRT_PI * GSL_res_array[legendre_array_index(k, 0, rot_inv)];
      }
      for (m = 1; m <= max_order; m++) {
        scaling_sin = M_2_SQRT_PI * M_SQRT2 * sin(-(m * phi));
        scaling_cos = M_2_SQRT_PI * M_SQRT2 * cos(m * phi);
        for (k = m; k <= max_order; k++) {
          (*sph)(i, Idxs2[k] - m) =
              scaling_sin * GSL_res_array[legendre_array_index(k, m, rot_inv)];
          (*sph)(i, Idxs2[k] + m) =
              scaling_cos * GSL_res_array[legendre_array_index(k, m, rot_inv)];
        }
      }
    }
  }

  return sph;
}

std::unique_ptr<Matrix<double>> spherical_bsplines1(
    const Matrix<double> &S, size_t n_basis,
    size_t degree,
    bool uniform_knot_angle_spacing) {
  auto basis = std::make_unique<Matrix<double>>(n_basis);
  std::vector<double> knots(n_basis + degree + 1);
  double s, s1, s2;
  std::vector<Matrix<double>> control(n_basis);
  std::vector<Matrix<double>> control_work(degree + 1);
  size_t interval;

  for (size_t i = 0; i <= degree; i++) {
    knots[i] = -1.0;
  }
  for (size_t i = degree + 1; i < n_basis; i++) {
    knots[i] = (double(i - degree) / double(n_basis - degree)) * 2.0 - 1.0;
    if (uniform_knot_angle_spacing) {
      knots[i] = sin(knots[i] * M_PI / 2.0);
    }
  }
  for (size_t i = n_basis; i <= n_basis + degree; i++) {
    knots[i] = 1.0;
  }

  for (size_t i = 0; i < n_basis; i++) {
    control[i] = Matrix<double>(n_basis);
    control[i](0, i) = 1.0;
  }

  FMLOG("degree\t" << degree << endl);
  FMLOG("n_basis\t" << n_basis << endl);
  FMLOG("n_basis+degree+1\t" << n_basis + degree + 1 << endl);

  for (size_t coord_idx = 0; coord_idx < S.rows(); coord_idx++) {
    s = S[coord_idx][0];

    FMLOG("step 1, coord_idx\t" << coord_idx << endl);
    interval = degree;
    while ((interval + 1 < n_basis) & (s >= knots[interval + 1]))
      interval++;

    FMLOG("step 2" << endl);
    for (size_t i = 0; i <= degree; i++)
      control_work[i] = control[i + interval - degree];

    FMLOG("step 3" << endl);
    for (size_t k = 1; k <= degree; k++) {
      for (size_t i = degree; i >= k; i--) {
        s1 = (knots[i + interval - k + 1] - s) /
          (knots[i + interval - k + 1] - knots[i + interval - degree]);
        s2 = 1.0 - s1;

        for (size_t j = 0; j < n_basis; j++)
          control_work[i](0, j) =
            (s1 * control_work[i - 1](0, j) + s2 * control_work[i](0, j));
      }
    }

    for (size_t j = 0; j < n_basis; j++) {
      (*basis)(coord_idx, j) = control_work[degree](0, j);
    }
  }

  return basis;
}

std::unique_ptr<Matrix<double>> spherical_bsplines(
    const Matrix3<double> &S, size_t n_basis,
    size_t degree,
    bool uniform_knot_angle_spacing) {
  auto basis = std::make_unique<Matrix<double>>(n_basis);
  std::vector<double> knots(n_basis + degree + 1);
  double s, s1, s2;
  std::vector<Matrix<double>> control(n_basis);
  std::vector<Matrix<double>> control_work(degree + 1);
  size_t interval;

  for (size_t i = 0; i <= degree; i++) {
    knots[i] = -1.0;
  }
  for (size_t i = degree + 1; i < n_basis; i++) {
    knots[i] = (double(i - degree) / double(n_basis - degree)) * 2.0 - 1.0;
    if (uniform_knot_angle_spacing) {
      knots[i] = sin(knots[i] * M_PI / 2.0);
    }
  }
  for (size_t i = n_basis; i <= n_basis + degree; i++) {
    knots[i] = 1.0;
  }

  for (size_t i = 0; i < n_basis; i++) {
    control[i] = Matrix<double>(n_basis);
    control[i](0, i) = 1.0;
  }

  FMLOG("degree\t" << degree << endl);
  FMLOG("n_basis\t" << n_basis << endl);
  FMLOG("n_basis+degree+1\t" << n_basis + degree + 1 << endl);

  for (size_t coord_idx = 0; coord_idx < S.rows(); coord_idx++) {
    s = S[coord_idx][2];

    FMLOG("step 1, coord_idx\t" << coord_idx << endl);
    interval = degree;
    while ((interval + 1 < n_basis) & (s >= knots[interval + 1]))
      interval++;

    FMLOG("step 2" << endl);
    for (size_t i = 0; i <= degree; i++)
      control_work[i] = control[i + interval - degree];

    FMLOG("step 3" << endl);
    for (size_t k = 1; k <= degree; k++) {
      for (size_t i = degree; i >= k; i--) {
        s1 = (knots[i + interval - k + 1] - s) /
          (knots[i + interval - k + 1] - knots[i + interval - degree]);
        s2 = 1.0 - s1;

        for (size_t j = 0; j < n_basis; j++)
          control_work[i](0, j) =
            (s1 * control_work[i - 1](0, j) + s2 * control_work[i](0, j));
      }
    }

    for (size_t j = 0; j < n_basis; j++) {
      (*basis)(coord_idx, j) = control_work[degree](0, j);
    }
  }

  return basis;
}

} /* namespace fmesh */
