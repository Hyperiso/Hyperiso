#ifndef __SCALAR_H__
#define __SCALAR_H__

#include <cmath>
#include "Utils.h"
#include "Logger.h"

template <class T>
std::enable_if_t<not std::numeric_limits<T>::is_integer, bool> fpeq(T, T, std::size_t n=10);

class scalar_t : public complex_t {
public:
    // scalar_t(double re = 0.0, double im = 0.0);
    constexpr scalar_t(double re = 0.0, double im = 0.0) : complex_t(re, im) {};
    scalar_t(complex_t z);
    
    scalar_t(const scalar_t& k);
    scalar_t& operator=(const scalar_t& k);
    scalar_t(scalar_t&& k) = default;
    scalar_t& operator=(scalar_t&& k) = default;

    operator double() const;

    scalar_t operator-() const;

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator+=(const T& rhs);

    friend scalar_t operator+(scalar_t lhs, const scalar_t& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator+(scalar_t lhs, const T& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator+(T lhs, const scalar_t& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator-=(const T& rhs);

    friend scalar_t operator-(scalar_t lhs, const scalar_t& rhs);
    
    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator-(scalar_t lhs, const T& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator-(T lhs, const scalar_t& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator*=(const T& rhs);

    friend scalar_t operator+(scalar_t lhs, const scalar_t& rhs);

    friend scalar_t operator+(double lhs, const scalar_t& rhs);

    friend scalar_t operator-(scalar_t lhs, const scalar_t& rhs);

    friend scalar_t operator-(double lhs, const scalar_t& rhs);

    friend scalar_t operator/(scalar_t lhs, const scalar_t& rhs);

    friend scalar_t operator/(double lhs, const scalar_t& rhs);

    friend scalar_t operator*(scalar_t lhs, const scalar_t& rhs);

    friend scalar_t operator*(double lhs, const scalar_t& rhs);

    friend scalar_t operator*(const scalar_t& lhs, double rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator*(scalar_t lhs, const T& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator*(T lhs, const scalar_t& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    scalar_t& operator/=(const T& rhs);

    friend scalar_t operator/(scalar_t lhs, const scalar_t& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator/(scalar_t lhs, const T& rhs);

    template<typename T>
    requires std::convertible_to<T, std::complex<double>>
    friend scalar_t operator/(T lhs, const scalar_t& rhs);

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

// TODO : warnings
// Define overload for pow(scalar_t, scalar_t)
scalar_t pow(const scalar_t& base, const scalar_t& exp);

// Define overload for pow(scalar_t, double)
scalar_t pow(const scalar_t& base, double exp);

template <typename T>
requires std::is_integral_v<T>
scalar_t pow(const scalar_t& base, T exp);

double real(const scalar_t& z);
double imag(const scalar_t& z);

// Ensure ADL works by defining functions in the same namespace
using std::sqrt;
using std::pow;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::exp;
using std::log;
using std::sinh;
using std::cosh;
using std::tanh;
using std::abs;
using std::arg;
using std::norm;
using std::real;
using std::imag;

#include "scalar.tpp"

#endif // __SCALAR_H__
