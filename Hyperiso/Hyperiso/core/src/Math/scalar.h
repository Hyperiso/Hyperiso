#ifndef __SCALAR_H__
#define __SCALAR_H__

#include <cmath>
#include "Utils.h"
#include "Logger.h"

template <class T>
std::enable_if_t<not std::numeric_limits<T>::is_integer, bool> fpeq(T, T, std::size_t n=10);

class scalar_t : public complex_t {
public:
    scalar_t(double re = 0.0, double im = 0.0) : complex_t(re, im) {};
    scalar_t(complex_t z) : complex_t(z) {};

    operator double() const {
        if (!fpeq(this->imag(), 0.)) {
            LOG_WARN("Casting complex values to double discards imaginary part.");
        }
        return this->real();
    };

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator+=(const T& rhs) {
        *this = static_cast<complex_t>(*this) + static_cast<complex_t>(rhs);
        return *this;
    }

    friend scalar_t operator+(scalar_t lhs, const scalar_t& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator+(scalar_t lhs, const T& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator+(T lhs, const scalar_t& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator-=(const T& rhs) {
        *this = static_cast<complex_t>(*this) - static_cast<complex_t>(rhs);
        return *this;
    }

    friend scalar_t operator-(scalar_t lhs, const scalar_t& rhs) {
        lhs -= rhs;
        return lhs;
    }
    
    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator-(scalar_t lhs, const T& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator-(T lhs, const scalar_t& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator*=(const T& rhs) {
        *this = static_cast<complex_t>(*this) * static_cast<complex_t>(rhs);
        return *this;
    }

    friend scalar_t operator*(scalar_t lhs, const scalar_t& rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator*(scalar_t lhs, const T& rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator*(T lhs, const scalar_t& rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator/=(const T& rhs) {
        *this = static_cast<complex_t>(*this) / static_cast<complex_t>(rhs);
        return *this;
    }

    friend scalar_t operator/(scalar_t lhs, const scalar_t& rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator/(scalar_t lhs, const T& rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator/(T lhs, const scalar_t& rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t pow(const scalar_t& base, T exponent) {
        return std::pow(static_cast<complex_t>(base), exponent);
    }
};
    
#define DEFINE_MATH_FUNCTION(func)                          \
inline scalar_t func(const scalar_t& x) {                   \
    return scalar_t(std::func(static_cast<complex_t>(x)));  \
}

DEFINE_MATH_FUNCTION(sqrt)
DEFINE_MATH_FUNCTION(sin)
DEFINE_MATH_FUNCTION(cos)
DEFINE_MATH_FUNCTION(tan)
DEFINE_MATH_FUNCTION(asin)
DEFINE_MATH_FUNCTION(acos)
DEFINE_MATH_FUNCTION(atan)
DEFINE_MATH_FUNCTION(exp)
DEFINE_MATH_FUNCTION(log)
DEFINE_MATH_FUNCTION(sinh)
DEFINE_MATH_FUNCTION(cosh)
DEFINE_MATH_FUNCTION(tanh)
DEFINE_MATH_FUNCTION(abs)
DEFINE_MATH_FUNCTION(arg)
DEFINE_MATH_FUNCTION(norm)

// Ensure ADL works by defining functions in the same namespace
namespace std {
    using ::sqrt;
    using ::pow;
    using ::sin;
    using ::cos;
    using ::tan;
    using ::asin;
    using ::acos;
    using ::atan;
    using ::exp;
    using ::log;
    using ::log10;
    using ::sinh;
    using ::cosh;
    using ::tanh;
    using ::abs;
    using ::arg;
    using ::norm;
}

#endif // __SCALAR_H__
