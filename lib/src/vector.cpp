// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <spamtrix_vector.hpp>
#include <spamtrix_exception.hpp>

namespace SpaMtrix {

namespace {

idx requireNonEmptyVector(const Vector& vector, const char* functionName) {
    if (vector.getLength() == 0) {
        throw SpaMtrixException("Vector must not be empty in " + std::string(functionName) + " at " + std::string(ERROR_LOCATION));
    }
    return vector.getLength();
}

}

Vector::Vector(const idx length) {
    values = std::vector<real>(length, 0.0);
    if (values.size() != length) {
      throw SpaMtrixException("error in " + std::string(ERROR_LOCATION) + "could not allocate " + std::to_string(length) + "elements");
    }
}

Vector::~Vector() {
}

Vector::Vector(const Vector &v) {
    values = v.values;
    if (values.size() != v.getLength()) {
      throw SpaMtrixException("error in " + std::string(ERROR_LOCATION) + "could not allocate " + std::to_string(v.getLength()) + "elements");
    }
}

Vector::Vector(const real *val, const idx &length) {
    values = std::vector<real>(val, val + length);
}

Vector &Vector::operator=(const Vector &v) {
    if (&v == this)
        return *this;
    values = v.values;
    return *this;
}

Vector &Vector::operator=(const real &a) {
    idx len = this->getLength();
    for (idx i = 0; i < len; ++i)
        values[i] = a;
    return *this;
}

void Vector::resize(const idx length) {
    values.resize(length, 0.0);
}

Vector &Vector::operator+=(const Vector &v) {
#ifdef DEBUG
    assert(this->getLength() == v.getLength());
#endif
    idx len = this->getLength();
    for (idx i = 0; i < len; ++i)
        this->values[i] += v.values[i];
    return *this;
}

Vector &Vector::operator-=(const Vector &v) {
#ifdef DEBUG
    assert(this->getLength() == v.getLength());
#endif
    idx len = this->getLength();
    for (idx i = 0; i < len; ++i)
        this->values[i] -= v.values[i];
    return *this;
}

Vector &Vector::operator+=(const real a) {
    idx len = this->getLength();
    for (idx i = 0; i < len; ++i)
        this->values[i] += a;
    return *this;
}

Vector &Vector::operator-=(const real a) {
    idx len = this->getLength();
    for (idx i = 0; i < len; ++i)
        this->values[i] -= a;
    return *this;
}

Vector &Vector::operator *=(const real a) {
    idx len = this->getLength();
    for (idx i = 0; i < len; ++i)
        this->values[i] *= a;
    return *this;
}

const Vector Vector::operator+(const Vector &rhs) const {
#ifdef DEBUG
    assert(this->getLength() == rhs.getLength());
#endif
    return Vector(*this) += rhs;
}

const Vector Vector::operator-(const Vector &rhs) const {
#ifdef DEBUG
    assert(this->getLength() == rhs.getLength());
#endif
    return Vector(*this) -= rhs;
}

const Vector Vector::operator*(const real a) const {
    return Vector(*this) *= a;
}

real Vector::getNorm() const {
    real sum(0);
    for (idx i = 0; i < this->getLength(); i++)
        sum += values[i] * values[i];
    return sqrt(sum);
}

real Vector::sum() const {
    real result(0);
    for (idx i = 0; i < this->getLength(); ++i) {
        result += values[i];
    }
    return result;
}

real Vector::max() const {
    return values[argMax()];
}

real Vector::min() const {
    return values[argMin()];
}

real Vector::absMax() const {
    return values[argAbsMax()];
}

real Vector::absMin() const {
    return values[argAbsMin()];
}

idx Vector::argMax() const {
    idx len = requireNonEmptyVector(*this, "argMax");
    idx bestIndex = 0;
    real bestValue = values[0];
    for (idx i = 1; i < len; ++i) {
        if (values[i] > bestValue) {
            bestValue = values[i];
            bestIndex = i;
        }
    }
    return bestIndex;
}

idx Vector::argMin() const {
    idx len = requireNonEmptyVector(*this, "argMin");
    idx bestIndex = 0;
    real bestValue = values[0];
    for (idx i = 1; i < len; ++i) {
        if (values[i] < bestValue) {
            bestValue = values[i];
            bestIndex = i;
        }
    }
    return bestIndex;
}

idx Vector::argAbsMax() const {
    idx len = requireNonEmptyVector(*this, "argAbsMax");
    idx bestIndex = 0;
    real bestValue = fabs(values[0]);
    for (idx i = 1; i < len; ++i) {
        real candidate = fabs(values[i]);
        if (candidate > bestValue) {
            bestValue = candidate;
            bestIndex = i;
        }
    }
    return bestIndex;
}

idx Vector::argAbsMin() const {
    idx len = requireNonEmptyVector(*this, "argAbsMin");
    idx bestIndex = 0;
    real bestValue = fabs(values[0]);
    for (idx i = 1; i < len; ++i) {
        real candidate = fabs(values[i]);
        if (candidate < bestValue) {
            bestValue = candidate;
            bestIndex = i;
        }
    }
    return bestIndex;
}

real Vector::maxAbsDiff(const Vector& other) const {
#ifdef DEBUG
    assert(this->getLength() == other.getLength());
#endif
    Vector diff(*this);
    if (diff.getLength() == 0) {
        return 0.0;
    }
    diff -= other;
    for (idx i = 0; i < diff.getLength(); ++i) {
        diff[i] = fabs(diff[i]);
    }
    return diff.max();
}

void Vector::normalise() {
    real norm = this->getNorm();
    if (norm == 0.0) // do nothing is zero-vector
        return;
    real k = 1.0 / norm;
    for (idx i = 0; i < this->getLength(); i++)
        values[i] *= k;
}
// Debug output.
void Vector::print(const char *name) const {
    std::cout << "Vector length " << this->getLength() << std::endl;
    for (idx i = 0; i < getLength(); i++) {
        if (name) {
            std::cout << name;
        } else {
            std::cout << "Vector";
        }
        std::cout << "[" << i << "] = " << values[i] << std::endl;
    }
}
} // end namespace SpaMtrix
