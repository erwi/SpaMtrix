// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>
#include <cmath>

#include "spamtrix_blas.hpp"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_tdmatrix.hpp"
#include "spamtrix_vector.hpp"

namespace {
using namespace SpaMtrix;

void requireVectorEquals(const Vector &actual, std::initializer_list<real> expected) {
  REQUIRE(actual.getLength() == expected.size());
  idx i = 0;
  for (real value : expected) {
    REQUIRE(actual[i] == Approx(value));
    ++i;
  }
}
}

TEST_CASE("BLAS level 1 operations behave correctly", "[BLAS]") {
  using namespace SpaMtrix;

  Vector v1(3);
  v1[0] = 1.0;
  v1[1] = -2.0;
  v1[2] = 3.0;

  Vector v2(3);
  v2[0] = 4.0;
  v2[1] = 5.0;
  v2[2] = 6.0;

  REQUIRE(dot(v1, v2) == Approx(12.0));
  REQUIRE(norm(v1) == Approx(std::sqrt(14.0)));

  Vector scaled = v1;
  scale(2.0, scaled);
  requireVectorEquals(scaled, {2.0, -4.0, 6.0});

  Vector y(3);
  y[0] = 1.0;
  y[1] = 1.0;
  y[2] = 1.0;
  axpy(2.0, v1, y);
  requireVectorEquals(y, {3.0, -3.0, 7.0});

  Vector z(3);
  z[0] = 2.0;
  z[1] = 4.0;
  z[2] = 6.0;
  aypx(0.5, z, v1);
  requireVectorEquals(z, {2.0, 0.0, 6.0});

  Vector magnitudes = abs(v1);
  requireVectorEquals(magnitudes, {1.0, 2.0, 3.0});
}

TEST_CASE("BLAS matrix-vector multiply and residual norms behave correctly", "[BLAS]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix I = mm.getIRCMatrix();

  Vector x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;

  Vector b(3);
  multiply(I, x, b);
  requireVectorEquals(b, {1.0, 2.0, 3.0});
  REQUIRE(multiply_dot(I, x, b) == Approx(14.0));

  REQUIRE(errorNorm2(I, x, b) < 1e-12);

  TDMatrix T(3);
  T.sparse_set(0, 0, 2.0);
  T.sparse_set(0, 1, -1.0);
  T.sparse_set(1, 0, -1.0);
  T.sparse_set(1, 1, 2.0);
  T.sparse_set(1, 2, -1.0);
  T.sparse_set(2, 1, -1.0);
  T.sparse_set(2, 2, 2.0);

  multiply(T, x, b);
  requireVectorEquals(b, {0.0, 0.0, 4.0});
  REQUIRE(errorNorm2(T, x, b) < 1e-12);
}