// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>

#include "spamtrix_densematrix.hpp"

namespace {
using namespace SpaMtrix;

void requireDenseEquals(const DenseMatrix &matrix, idx rows, idx cols, real expected) {
  REQUIRE(matrix.getNumRows() == rows);
  REQUIRE(matrix.getNumCols() == cols);
  for (idx c = 0; c < cols; ++c) {
    for (idx r = 0; r < rows; ++r) {
      REQUIRE(matrix(r, c) == Approx(expected));
    }
  }
}
}

TEST_CASE("DenseMatrix construction initializes all entries to zero", "[DenseMatrix]") {
  using namespace SpaMtrix;

  DenseMatrix matrix(3, 4);
  requireDenseEquals(matrix, 3, 4, 0.0);
}

TEST_CASE("DenseMatrix operator() supports read and write access", "[DenseMatrix]") {
  using namespace SpaMtrix;

  DenseMatrix matrix(2, 3);
  matrix(0, 0) = 1.0;
  matrix(1, 0) = 2.0;
  matrix(0, 1) = 3.0;
  matrix(1, 2) = 4.0;

  REQUIRE(matrix(0, 0) == Approx(1.0));
  REQUIRE(matrix(1, 0) == Approx(2.0));
  REQUIRE(matrix(0, 1) == Approx(3.0));
  REQUIRE(matrix(1, 2) == Approx(4.0));
}

TEST_CASE("DenseMatrix setAllValuesTo overwrites every stored value", "[DenseMatrix]") {
  using namespace SpaMtrix;

  DenseMatrix matrix(2, 2);
  matrix(0, 0) = 1.0;
  matrix(1, 1) = 2.0;

  matrix.setAllValuesTo(5.5);
  requireDenseEquals(matrix, 2, 2, 5.5);
}

TEST_CASE("DenseMatrix copy construction produces an independent copy", "[DenseMatrix]") {
  using namespace SpaMtrix;

  DenseMatrix original(2, 2);
  original(0, 0) = 1.0;
  original(1, 1) = 2.0;

  DenseMatrix copy(original);
  original(0, 0) = 9.0;
  original(1, 1) = 8.0;

  REQUIRE(copy(0, 0) == Approx(1.0));
  REQUIRE(copy(1, 1) == Approx(2.0));
}

TEST_CASE("DenseMatrix assignment produces an independent copy", "[DenseMatrix]") {
  using namespace SpaMtrix;

  DenseMatrix source(2, 2);
  source(0, 1) = 3.5;

  DenseMatrix target(1, 1);
  target = source;

  REQUIRE(target.getNumRows() == 2);
  REQUIRE(target.getNumCols() == 2);
  REQUIRE(target(0, 1) == Approx(3.5));

  source(0, 1) = 7.5;
  REQUIRE(target(0, 1) == Approx(3.5));
}