#include "scalar.h"

// scalar_t::scalar_t(double re, double im) : complex_t(re, im) {}
// constexpr scalar_t::scalar_t(double re, double im) 

scalar_t::scalar_t(complex_t z) : complex_t(z) {}

scalar_t::scalar_t(const scalar_t& k) : complex_t(static_cast<complex_t>(k)) {}

scalar_t& scalar_t::operator=(const scalar_t& k) {
    complex_t::operator=(static_cast<const complex_t&>(k));
    return *this;
}

scalar_t::operator double() const {
    if (!fpeq(this->imag(), 0.)) {
        LOG_WARN("Casting complex values to double discards imaginary part: (", this->real(), ",", this->imag(), ")");
    }
    return this->real();
};

scalar_t scalar_t::operator-() const {
    return -static_cast<complex_t>(*this);
}

scalar_t operator+(double lhs, const scalar_t& rhs) {
    return scalar_t(lhs) + rhs;
}

scalar_t operator+(const scalar_t& lhs, double rhs) {
    return lhs + scalar_t(rhs);
}

scalar_t operator+(scalar_t lhs, const scalar_t& rhs) {
    lhs += rhs;
    return lhs;
}

scalar_t operator-(scalar_t lhs, const scalar_t& rhs) {
    lhs -= rhs;
    return lhs;
}

scalar_t operator-(double lhs, const scalar_t& rhs) {
    return scalar_t(lhs) - rhs;
}

scalar_t operator-(const scalar_t& lhs, double rhs) {
    return lhs - scalar_t(rhs);
}

scalar_t operator*(scalar_t lhs, const scalar_t& rhs) {
    lhs *= rhs;
    return lhs;
}

scalar_t operator*(double lhs, const scalar_t& rhs) {
    return scalar_t(lhs) * rhs;
}

scalar_t operator*(const scalar_t& lhs, double rhs) {
    return lhs * scalar_t(rhs);
}

scalar_t operator/(scalar_t lhs, const scalar_t& rhs) {
    lhs /= rhs;
    return lhs;
}

scalar_t operator/(double lhs, const scalar_t& rhs) {
    return scalar_t(lhs) / rhs;
}

scalar_t operator/(const scalar_t& lhs, double rhs) {
    return lhs / scalar_t(rhs);
}

scalar_t pow(const scalar_t& base, const scalar_t& exp) {
    return scalar_t(std::pow(static_cast<complex_t>(base), static_cast<complex_t>(exp)));
}

scalar_t pow(const scalar_t& base, double exp) {
    return scalar_t(std::pow(static_cast<complex_t>(base), exp));
}

double real(const scalar_t& z) { 
    return z.real(); 
}

double imag(const scalar_t& z) { 
    return z.imag(); 
}