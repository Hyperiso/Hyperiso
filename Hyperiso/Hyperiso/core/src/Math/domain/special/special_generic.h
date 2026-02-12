#ifndef __SPECIAL_GENERIC_H__
#define __SPECIAL_GENERIC_H__

#include <cmath>
#include "../linalg/scalar.h"

/**
 * @brief Computes the exponential integral function Ei(x).
 * @param x Input value.
 * @return Result of the exponential integral.
 */
double Ei(double x);

/**
 * @brief Computes the dilogarithm function Li2(x).
 * @param x Input value.
 * @return Result of the dilogarithm.
 */
double Li2(double x);

/**
 * @brief Computes the complex dilogarithm function.
 * @param x Complex input value.
 * @return Complex dilogarithm result.
 */
scalar_t CLi2(scalar_t x);

/**
 * @brief Computes the complex trilogarithm function.
 * @param x Complex input value.
 * @return Complex trilogarithm result.
 */
scalar_t CLi3(scalar_t x);

/**
 * @brief Computes the complex quadrilogarithm function.
 * @param x Complex input value.
 * @return Complex quadrilogarithm result.
 */
scalar_t CLi4(scalar_t x);

/**
 * @brief Computes the trilogarithm function Li3(x).
 * @param x Input value.
 * @return Result of the trilogarithm.
 */
double Li3(double x);

/**
 * @brief Computes the Clausen function Cl2(x).
 * @param x Input value.
 * @return Result of the Clausen function.
 */
double Cl2(double x);

/**
 * @brief Computes the Clausen function Cl3(x).
 * @param x Input value.
 * @return Result of the Clausen function.
 */
double Cl3(double x);

// Bessel

double I0(double x);
double I1(double x);
double K0(double x);
double K1(double x);
double K2(double x);
double K3(double x);
double K4(double x);
double Lbessel(double x);
double Mbessel(double x);
double Nbessel(double x);
double K0exp(double x, double z);
double K1exp(double x, double z);
double K2exp(double x, double z);

/**
 * @brief Kronecker delta function.
 * @param x First integer.
 * @param y Second integer.
 * @return 1 if x == y, else 0.
 */
double kron(int x, int y);

/**
 * @brief Returns the sign of a value.
 * @tparam T Type (must be comparable).
 * @param val Value to check.
 * @return -1, 0, or 1 depending on sign of val.
 */
template <typename T> int sgn(T val) {
    return (T(0) < val) -(val<T(0));
}

/**
 * @brief Compares two floating point numbers with a given precision.
 * @tparam T Floating point type.
 * @param x First number.
 * @param y Second number.
 * @param n Number of machine epsilons allowed.
 * @return True if numbers are approximately equal.
 */
template <class T>
std::enable_if_t<not std::numeric_limits<T>::is_integer, bool>
fpeq(T x, T y, std::size_t n) {
    const T m = std::min(std::fabs(x), std::fabs(y));
    const int exp = m < std::numeric_limits<T>::min()
                  ? std::numeric_limits<T>::min_exponent - 1
                  : std::ilogb(m);
    return std::fabs(x - y) <= n * std::ldexp(std::numeric_limits<T>::epsilon(), exp);
}

#endif // __SPECIAL_GENERIC_H__
