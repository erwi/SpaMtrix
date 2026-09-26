// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <catch.h>
#include <spamtrix_vector.hpp>
#include <spamtrix_exception.hpp>

TEST_CASE("Vector construction initializes all entries to zero", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(4);

    REQUIRE(v.getLength() == 4);
    for (idx i = 0; i < v.getLength(); ++i) {
        REQUIRE(v[i] == 0.0);
    }
}

TEST_CASE("Vector copy construction creates an independent copy", "[Vector]") {
    using namespace SpaMtrix;

    Vector original(3);
    original[0] = 1.0;
    original[1] = 2.0;
    original[2] = 3.0;

    Vector copy(original);
    original[1] = 42.0;

    REQUIRE(copy.getLength() == 3);
    REQUIRE(copy[0] == 1.0);
    REQUIRE(copy[1] == 2.0);
    REQUIRE(copy[2] == 3.0);
}

TEST_CASE("Vector assignment creates an independent deep copy", "[Vector]") {
    using namespace SpaMtrix;

    Vector source(3);
    source[0] = 4.0;
    source[1] = 5.0;
    source[2] = 6.0;

    Vector target(3);
    target = source;
    source[2] = -1.0;

    REQUIRE(target.getLength() == 3);
    REQUIRE(target[0] == 4.0);
    REQUIRE(target[1] == 5.0);
    REQUIRE(target[2] == 6.0);
}

TEST_CASE("Vector scalar assignment sets every element", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(5);
    v = 2.5;

    for (idx i = 0; i < v.getLength(); ++i) {
        REQUIRE(v[i] == 2.5);
    }
}

TEST_CASE("Vector element access supports read and write", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(2);
    v[0] = 7.0;
    v[1] = 11.0;

    REQUIRE(v[0] == 7.0);
    REQUIRE(v[1] == 11.0);
    REQUIRE(v(0) == 7.0);
    REQUIRE(v(1) == 11.0);
}

TEST_CASE("Vector compound arithmetic with scalars and vectors works", "[Vector]") {
    using namespace SpaMtrix;

    Vector v1(3);
    v1[0] = 1.0;
    v1[1] = 2.0;
    v1[2] = 3.0;

    Vector v2(3);
    v2[0] = 4.0;
    v2[1] = 5.0;
    v2[2] = 6.0;

    v1 += v2;
    REQUIRE(v1[0] == 5.0);
    REQUIRE(v1[1] == 7.0);
    REQUIRE(v1[2] == 9.0);

    v1 -= v2;
    REQUIRE(v1[0] == 1.0);
    REQUIRE(v1[1] == 2.0);
    REQUIRE(v1[2] == 3.0);

    v1 += 1.0;
    REQUIRE(v1[0] == 2.0);
    REQUIRE(v1[1] == 3.0);
    REQUIRE(v1[2] == 4.0);

    v1 -= 2.0;
    REQUIRE(v1[0] == 0.0);
    REQUIRE(v1[1] == 1.0);
    REQUIRE(v1[2] == 2.0);

    v1 *= 3.0;
    REQUIRE(v1[0] == 0.0);
    REQUIRE(v1[1] == 3.0);
    REQUIRE(v1[2] == 6.0);
}

TEST_CASE("Vector non-mutating arithmetic returns new vectors", "[Vector]") {
    using namespace SpaMtrix;

    Vector v1(3);
    v1[0] = 1.0;
    v1[1] = 2.0;
    v1[2] = 3.0;

    Vector v2(3);
    v2[0] = 4.0;
    v2[1] = 5.0;
    v2[2] = 6.0;

    Vector sum = v1 + v2;
    Vector diff = v2 - v1;
    Vector scaled_right = v1 * 2.0;
    Vector scaled_left = 2.0 * v1;

    REQUIRE(v1[0] == 1.0);
    REQUIRE(v1[1] == 2.0);
    REQUIRE(v1[2] == 3.0);

    REQUIRE(sum[0] == 5.0);
    REQUIRE(sum[1] == 7.0);
    REQUIRE(sum[2] == 9.0);

    REQUIRE(diff[0] == 3.0);
    REQUIRE(diff[1] == 3.0);
    REQUIRE(diff[2] == 3.0);

    REQUIRE(scaled_right[0] == 2.0);
    REQUIRE(scaled_right[1] == 4.0);
    REQUIRE(scaled_right[2] == 6.0);

    REQUIRE(scaled_left[0] == 2.0);
    REQUIRE(scaled_left[1] == 4.0);
    REQUIRE(scaled_left[2] == 6.0);
}

TEST_CASE("Vector setAllValuesTo, norm, and normalise behave correctly", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(2);
    v[0] = 3.0;
    v[1] = 4.0;

    REQUIRE(v.getNorm() == Approx(5.0));

    v.normalise();
    REQUIRE(v.getNorm() == Approx(1.0));
    REQUIRE(v[0] == Approx(0.6));
    REQUIRE(v[1] == Approx(0.8));

    v.setAllValuesTo(7.0);
    REQUIRE(v[0] == 7.0);
    REQUIRE(v[1] == 7.0);
}

TEST_CASE("Vector reduction helpers return expected values and indices", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(5);
    v[0] = -2.0;
    v[1] = 4.0;
    v[2] = -7.5;
    v[3] = 0.5;
    v[4] = 4.0;

    REQUIRE(v.sum() == Approx(-1.0));
    REQUIRE(v.max() == Approx(4.0));
    REQUIRE(v.min() == Approx(-7.5));
    REQUIRE(v.absMax() == Approx(-7.5));
    REQUIRE(v.absMin() == Approx(0.5));

    REQUIRE(v.argMax() == 1);
    REQUIRE(v.argMin() == 2);
    REQUIRE(v.argAbsMax() == 2);
    REQUIRE(v.argAbsMin() == 3);
}

TEST_CASE("Vector reduction helpers prefer the first index on ties", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(4);
    v[0] = -3.0;
    v[1] = 3.0;
    v[2] = -3.0;
    v[3] = 3.0;

    REQUIRE(v.argMax() == 1);
    REQUIRE(v.argMin() == 0);
    REQUIRE(v.argAbsMax() == 0);
}

TEST_CASE("Vector maxAbsDiff returns the largest element-wise difference", "[Vector]") {
    using namespace SpaMtrix;

    Vector lhs(4);
    lhs[0] = 1.0;
    lhs[1] = -2.0;
    lhs[2] = 3.5;
    lhs[3] = 0.0;

    Vector rhs(4);
    rhs[0] = -1.0;
    rhs[1] = -5.5;
    rhs[2] = 3.0;
    rhs[3] = 2.0;

    REQUIRE(lhs.maxAbsDiff(rhs) == Approx(3.5));
    REQUIRE(rhs.maxAbsDiff(lhs) == Approx(3.5));
}

TEST_CASE("Vector reductions on an empty vector throw", "[Vector]") {
    using namespace SpaMtrix;

    Vector empty;

    REQUIRE_THROWS_AS(empty.max(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.min(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.absMax(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.absMin(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.argMax(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.argMin(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.argAbsMax(), SpaMtrixException);
    REQUIRE_THROWS_AS(empty.argAbsMin(), SpaMtrixException);
    REQUIRE(empty.sum() == Approx(0.0));
    REQUIRE(empty.maxAbsDiff(empty) == Approx(0.0));
}

TEST_CASE("Vector resize preserves existing values and zero-fills new elements", "[Vector]") {
    using namespace SpaMtrix;

    Vector v(3);
    v[0] = 1.0;
    v[1] = 2.0;
    v[2] = 3.0;

    v.resize(5);
    REQUIRE(v.getLength() == 5);
    REQUIRE(v[0] == 1.0);
    REQUIRE(v[1] == 2.0);
    REQUIRE(v[2] == 3.0);
    REQUIRE(v[3] == 0.0);
    REQUIRE(v[4] == 0.0);

    v[3] = 9.0;
    v[4] = 10.0;
    v.resize(2);
    REQUIRE(v.getLength() == 2);
    REQUIRE(v[0] == 1.0);
    REQUIRE(v[1] == 2.0);
}