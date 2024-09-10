#include <catch.h>
#include "spamtrix_matrixmaker.hpp"

TEST_CASE("Add sparse matrix to other - identical sparsity patterns") {
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

TEST_CASE("Add sparse matrix to other - different sparsity patterns") {
  using namespace SpaMtrix;
  // ARRANGE
  MatrixMaker mm1(3, 3);
  mm1.identity();
  auto I1 = mm1.getIRCMatrix();

  MatrixMaker mm2(2, 2);
  mm2.addNonZero(1, 1, 1);
  auto M = mm2.getIRCMatrix();

  // ACT
  I1.add(M);

  // ASSERT
  REQUIRE(I1.getValue(0, 0) == 1);
  REQUIRE(I1.getValue(1, 1) == 2);
  REQUIRE(I1.getValue(2, 2) == 1);
}



