// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>

#include "spamtrix_fleximatrix.hpp"
#include "spamtrix_matrixmaker.hpp"

namespace {
using namespace SpaMtrix;

void requireFlexiValue(const FlexiMatrix &matrix, idx row, idx col, real expected) {
  REQUIRE(matrix.getValue(row, col) == Approx(expected));
}
}

TEST_CASE("FlexiMatrix addNonZero stores values and calcNumNonZeros counts them", "[FlexiMatrix]") {
  using namespace SpaMtrix;

  FlexiMatrix matrix;
  matrix.addNonZero(0, 0, 1.0);
  matrix.addNonZero(0, 2, 3.0);
  matrix.addNonZero(1, 1, 2.0);

  REQUIRE(matrix.getNumRows() == 2);
  REQUIRE(matrix.getNumCols() == 3);
  REQUIRE(matrix.calcNumNonZeros() == 3);
  requireFlexiValue(matrix, 0, 0, 1.0);
  requireFlexiValue(matrix, 0, 2, 3.0);
  requireFlexiValue(matrix, 1, 1, 2.0);
  REQUIRE(matrix.getValue(0, 1) == Approx(0.0));
}

TEST_CASE("FlexiMatrix isNonZero returns a pointer to stored data", "[FlexiMatrix]") {
  using namespace SpaMtrix;

  FlexiMatrix matrix;
  matrix.addNonZero(1, 1, 4.0);

  real *valuePtr = nullptr;
  REQUIRE(matrix.isNonZero(1, 1, valuePtr));
  REQUIRE(valuePtr != nullptr);
  REQUIRE(*valuePtr == Approx(4.0));

  valuePtr = nullptr;
  REQUIRE_FALSE(matrix.isNonZero(0, 0, valuePtr));
  REQUIRE(valuePtr == nullptr);
}

TEST_CASE("FlexiMatrix setValue updates existing storage and inserts new storage", "[FlexiMatrix]") {
  using namespace SpaMtrix;

  FlexiMatrix matrix;
  matrix.addNonZero(0, 0, 1.0);

  matrix.setValue(0, 0, 2.0);
  requireFlexiValue(matrix, 0, 0, 2.0);

  matrix.setValue(0, 1, 3.0);
  requireFlexiValue(matrix, 0, 1, 3.0);
  REQUIRE(matrix.calcNumNonZeros() == 2);
}

TEST_CASE("FlexiMatrix construction from IRCMatrix preserves non-zeros", "[FlexiMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix source = mm.getIRCMatrix();

  FlexiMatrix matrix(source);

  REQUIRE(matrix.getNumRows() == 3);
  REQUIRE(matrix.getNumCols() == 3);
  REQUIRE(matrix.calcNumNonZeros() == 3);
  requireFlexiValue(matrix, 0, 0, 1.0);
  requireFlexiValue(matrix, 1, 1, 1.0);
  requireFlexiValue(matrix, 2, 2, 1.0);
}