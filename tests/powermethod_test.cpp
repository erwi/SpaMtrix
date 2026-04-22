// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>

#include <cmath>

#include "spamtrix_ircmatrix.hpp"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_powermethod.hpp"
#include "spamtrix_vector.hpp"

namespace {
using namespace SpaMtrix;
}

TEST_CASE("PowerMethod finds the dominant eigenpair of the identity matrix", "[PowerMethod]") {
  using namespace SpaMtrix;

  MatrixMaker mm(4, 4);
  mm.identity();
  IRCMatrix A = mm.getIRCMatrix();

  Vector eigenVector(4);
  real eigenValue = 0.0;
  real toler = 1e-12;

  idx iterations = powerMethod(A, eigenValue, eigenVector, toler, 10);

  REQUIRE(iterations > 0);
  REQUIRE(eigenValue == Approx(1.0));
  REQUIRE(eigenVector[0] == Approx(1.0));
  REQUIRE(eigenVector[1] == Approx(0.0));
  REQUIRE(eigenVector[2] == Approx(0.0));
  REQUIRE(eigenVector[3] == Approx(0.0));
}

TEST_CASE("PowerMethod converges to the largest diagonal entry", "[PowerMethod]") {
  using namespace SpaMtrix;

  MatrixMaker mm(4, 4);
  mm.addNonZero(0, 0, 10.0);
  mm.addNonZero(1, 1, 5.0);
  mm.addNonZero(2, 2, 2.0);
  mm.addNonZero(3, 3, 1.0);
  IRCMatrix A = mm.getIRCMatrix();

  Vector eigenVector(4);
  eigenVector[0] = 1.0;
  eigenVector[1] = 1.0;
  eigenVector[2] = 1.0;
  eigenVector[3] = 1.0;

  real eigenValue = 0.0;
  real toler = 1e-10;

  idx iterations = powerMethod(A, eigenValue, eigenVector, toler, 25);

  REQUIRE(iterations > 0);
  REQUIRE(iterations <= 25);
  REQUIRE(eigenValue == Approx(10.0));
  REQUIRE(std::abs(eigenVector[0]) > 0.99);
  REQUIRE(std::abs(eigenVector[1]) < 0.1);
  REQUIRE(std::abs(eigenVector[2]) < 0.1);
  REQUIRE(std::abs(eigenVector[3]) < 0.1);
}

TEST_CASE("PowerMethod respects the maxIter limit", "[PowerMethod]") {
  using namespace SpaMtrix;

  MatrixMaker mm(4, 4);
  mm.addNonZero(0, 0, 10.0);
  mm.addNonZero(1, 1, 5.0);
  mm.addNonZero(2, 2, 2.0);
  mm.addNonZero(3, 3, 1.0);
  IRCMatrix A = mm.getIRCMatrix();

  Vector eigenVector(4);
  eigenVector[0] = 1.0;
  eigenVector[1] = 1.0;
  eigenVector[2] = 1.0;
  eigenVector[3] = 1.0;

  real eigenValue = 0.0;
  real toler = 1e-14;

  idx iterations = powerMethod(A, eigenValue, eigenVector, toler, 1);

  REQUIRE(iterations == 1);
}