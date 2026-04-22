// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>
#include <utility>

#include "spamtrix_blas.hpp"
#include "spamtrix_exception.hpp"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_vector.hpp"

namespace {
using namespace SpaMtrix;

void requireMatrixValue(const IRCMatrix &matrix, idx row, idx col, real expected) {
  REQUIRE(matrix.getValue(row, col) == Approx(expected));
}
}

TEST_CASE("Add sparse matrix to other - identical sparsity patterns", "[IRCMatrix]") {
  using namespace SpaMtrix;
  MatrixMaker mm(3, 3);
  mm.identity();
  auto I1 = mm.getIRCMatrix();
  auto I2 = mm.getIRCMatrix();


  SECTION("Scale by 1") {
    I2.add(I1);

    REQUIRE(I2.getValue(0, 0) == 2);
    REQUIRE(I2.getValue(1, 1) == 2);
    REQUIRE(I2.getValue(2, 2) == 2);
  }

  SECTION("Scale by 2") {
    I2.add(I1, 2);

    REQUIRE(I2.getValue(0, 0) == 3);
    REQUIRE(I2.getValue(1, 1) == 3);
    REQUIRE(I2.getValue(2, 2) == 3);
  }
}

TEST_CASE("IRCMatrix vector multiplication matches identity action", "[IRCMatrix][multiply]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix A = mm.getIRCMatrix();

  Vector x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;

  Vector y = A * x;

  REQUIRE(y.getLength() == 3);
  REQUIRE(y[0] == Approx(1.0));
  REQUIRE(y[1] == Approx(2.0));
  REQUIRE(y[2] == Approx(3.0));

  Vector b(3);
  multiply(A, x, b);
  REQUIRE(b[0] == Approx(1.0));
  REQUIRE(b[1] == Approx(2.0));
  REQUIRE(b[2] == Approx(3.0));
}

TEST_CASE("IRCMatrix scalar multiplication updates values in-place and by copy", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix A = mm.getIRCMatrix();

  IRCMatrix scaled = A * 4.0;
  requireMatrixValue(A, 0, 0, 1.0);
  requireMatrixValue(A, 1, 1, 1.0);
  requireMatrixValue(A, 2, 2, 1.0);
  requireMatrixValue(scaled, 0, 0, 4.0);
  requireMatrixValue(scaled, 1, 1, 4.0);
  requireMatrixValue(scaled, 2, 2, 4.0);

  A *= 2.0;
  requireMatrixValue(A, 0, 0, 2.0);
  requireMatrixValue(A, 1, 1, 2.0);
  requireMatrixValue(A, 2, 2, 2.0);
}

TEST_CASE("IRCMatrix sparse set, add, and get round-trip existing storage", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix A = mm.getIRCMatrix();

  A.sparse_set(0, 0, 5.0);
  REQUIRE(A.sparse_get(0, 0) == Approx(5.0));

  A.sparse_add(1, 1, 2.0);
  REQUIRE(A.sparse_get(1, 1) == Approx(3.0));

  REQUIRE_THROWS_AS(A.sparse_set(0, 1, 7.0), SpaMtrixException);
  REQUIRE_THROWS_AS(A.sparse_add(0, 1, 1.0), SpaMtrixException);
  REQUIRE_THROWS_AS(A.sparse_get(0, 1), SpaMtrixException);
}

TEST_CASE("IRCMatrix copy construction produces an independent deep copy", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix original = mm.getIRCMatrix();
  IRCMatrix copy(original);

  copy.sparse_set(0, 0, 9.0);

  requireMatrixValue(original, 0, 0, 1.0);
  requireMatrixValue(copy, 0, 0, 9.0);
}

TEST_CASE("IRCMatrix move construction transfers ownership and clears source", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix original = mm.getIRCMatrix();

  IRCMatrix moved(std::move(original));

  REQUIRE(original.getNumRows() == 0);
  REQUIRE(original.getNumCols() == 0);
  REQUIRE(original.getnnz() == 0);

  requireMatrixValue(moved, 0, 0, 1.0);
  requireMatrixValue(moved, 1, 1, 1.0);
  requireMatrixValue(moved, 2, 2, 1.0);
}

TEST_CASE("IRCMatrix move assignment transfers ownership and clears source", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix source = mm.getIRCMatrix();

  IRCMatrix target;
  target = std::move(source);

  REQUIRE(source.getNumRows() == 0);
  REQUIRE(source.getNumCols() == 0);
  REQUIRE(source.getnnz() == 0);

  requireMatrixValue(target, 0, 0, 1.0);
  requireMatrixValue(target, 1, 1, 1.0);
  requireMatrixValue(target, 2, 2, 1.0);
}

TEST_CASE("IRCMatrix scalar assignment sets all stored values", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix A = mm.getIRCMatrix();

  A = 2.5;

  requireMatrixValue(A, 0, 0, 2.5);
  requireMatrixValue(A, 1, 1, 2.5);
  requireMatrixValue(A, 2, 2, 2.5);
}

TEST_CASE("IRCMatrix getValuePtr and multiply helper remain consistent", "[IRCMatrix][multiply]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  IRCMatrix A = mm.getIRCMatrix();

  REQUIRE(A.getValuePtr(1, 1) != nullptr);
  REQUIRE(*A.getValuePtr(1, 1) == Approx(1.0));

  Vector x(3);
  x[0] = 4.0;
  x[1] = 5.0;
  x[2] = 6.0;
  Vector b(3);
  multiply(A, x, b);

  REQUIRE(b[0] == Approx(4.0));
  REQUIRE(b[1] == Approx(5.0));
  REQUIRE(b[2] == Approx(6.0));
}

TEST_CASE("Add sparse matrix to other - different sparsity patterns", "[IRCMatrix]") {
  using namespace SpaMtrix;
  // Arrange.
  MatrixMaker mm1(3, 3);
  mm1.identity();
  auto I1 = mm1.getIRCMatrix();

  MatrixMaker mm2(2, 2);
  mm2.addNonZero(1, 1, 1);
  auto M = mm2.getIRCMatrix();

  // Act.
  I1.add(M);

  // Assert.
  REQUIRE(I1.getValue(0, 0) == 1);
  REQUIRE(I1.getValue(1, 1) == 2);
  REQUIRE(I1.getValue(2, 2) == 1);
}

TEST_CASE("Access matrix values by pointer", "[IRCMatrix]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 3);
  mm.identity();
  auto I = mm.getIRCMatrix();

  REQUIRE(I.getValuePtr(0, 0) != nullptr);
  REQUIRE(*I.getValuePtr(1, 1) == 1.);

  // Modify the matrix value.
  *I.getValuePtr(2, 2) = 2;
  REQUIRE(*I.getValuePtr(2, 2) == 2);
}

