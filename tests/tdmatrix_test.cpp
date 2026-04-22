// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>

#include "spamtrix_blas.hpp"
#include "spamtrix_tdmatrix.hpp"
#include "spamtrix_vector.hpp"

namespace {
using namespace SpaMtrix;

void requireVectorEquals(const Vector &actual, const std::initializer_list<real> &expected) {
  REQUIRE(actual.getLength() == expected.size());
  idx i = 0;
  for (real value : expected) {
    REQUIRE(actual[i] == Approx(value));
    ++i;
  }
}
}

TEST_CASE("TDMatrix construction reports the expected dimensions and nnz", "[TDMatrix]") {
  using namespace SpaMtrix;

  TDMatrix A(4);

  REQUIRE(A.getNumRows() == 4);
  REQUIRE(A.getNumCols() == 4);
  REQUIRE(A.getnnz() == 10);
}

TEST_CASE("TDMatrix sparse set and get round-trip on all bands", "[TDMatrix]") {
  using namespace SpaMtrix;

  TDMatrix A(3);
  A.sparse_set(0, 0, 2.0);
  A.sparse_set(0, 1, -1.0);
  A.sparse_set(1, 0, -1.0);
  A.sparse_set(1, 1, 2.0);
  A.sparse_set(1, 2, -1.0);
  A.sparse_set(2, 1, -1.0);
  A.sparse_set(2, 2, 2.0);

  REQUIRE(A.sparse_get(0, 0) == Approx(2.0));
  REQUIRE(A.sparse_get(0, 1) == Approx(-1.0));
  REQUIRE(A.sparse_get(1, 0) == Approx(-1.0));
  REQUIRE(A.sparse_get(1, 1) == Approx(2.0));
  REQUIRE(A.sparse_get(1, 2) == Approx(-1.0));
  REQUIRE(A.sparse_get(2, 1) == Approx(-1.0));
  REQUIRE(A.sparse_get(2, 2) == Approx(2.0));
}

TEST_CASE("TDMatrix multiply matches a known tridiagonal stencil", "[TDMatrix][BLAS]") {
  using namespace SpaMtrix;

  TDMatrix A(3);
  A.sparse_set(0, 0, 2.0);
  A.sparse_set(0, 1, -1.0);
  A.sparse_set(1, 0, -1.0);
  A.sparse_set(1, 1, 2.0);
  A.sparse_set(1, 2, -1.0);
  A.sparse_set(2, 1, -1.0);
  A.sparse_set(2, 2, 2.0);

  Vector x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;

  Vector b(3);
  multiply(A, x, b);

  requireVectorEquals(b, {0.0, 0.0, 4.0});
}

TEST_CASE("TDMatrix solveAxb recovers an exact solution", "[TDMatrix][solver]") {
  using namespace SpaMtrix;

  TDMatrix A(3);
  A.sparse_set(0, 0, 2.0);
  A.sparse_set(0, 1, -1.0);
  A.sparse_set(1, 0, -1.0);
  A.sparse_set(1, 1, 2.0);
  A.sparse_set(1, 2, -1.0);
  A.sparse_set(2, 1, -1.0);
  A.sparse_set(2, 2, 2.0);

  Vector expected(3);
  expected[0] = 1.0;
  expected[1] = 2.0;
  expected[2] = 3.0;

  Vector b(3);
  multiply(A, expected, b);

  Vector x(3);
  A.solveAxb(x, b);

  requireVectorEquals(x, {1.0, 2.0, 3.0});
  REQUIRE(errorNorm2(A, x, b) < 1e-12);
}