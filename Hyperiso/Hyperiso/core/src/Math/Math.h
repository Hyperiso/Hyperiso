#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <type_traits>
#include "Utils.h"
#include "SparseMatrix.h"
#include "Logger.h"
#include "scalar.h"
#include "analysis.h"

/* Constants */

constexpr double PI =           3.1415926535897932;
constexpr double PI2 =          9.869604401089358;
constexpr double PI3 =          31.006276680299816;
constexpr double INV_PI =       0.318309886183791;
constexpr double INV_PI2 =      0.101321183642338;
constexpr double INV_PI3 =      0.032251534433199;
constexpr double E =            2.7182818284590452;
constexpr double ZETA3 =        1.2020569031595942;
constexpr double RT2 =          1.414213562373095;
constexpr double INV_RT2 =      0.707106781186547;
constexpr double ERF_INV_RT2 =  0.682689492137085;

constexpr double HBAR =     6.58211889e-25;  // GeV.s
constexpr double G_NEWT =   6.67428e-8;  // cm^3.g^-1.s^-2
constexpr double M_PLANCK = 1.2209102930946623e+19; // GeV
constexpr double M_NUCL =   0.939;  // Gev
constexpr double M_P =   0.9382720;  // Gev
constexpr double M_N =   0.9395654;  // Gev

constexpr double EPSILON = 1e-5;

constexpr complex_t I = complex_t(0, 1);

/* Conversion factors */

constexpr double GEV_TO_INV_PS = 1.519268e12;
constexpr double GEV_TO_INV_S = 1.519268e24;

/* Functions */

/**
 * @brief Kronecker delta function.
 * @param x First integer.
 * @param y Second integer.
 * @return 1 if x == y, else 0.
 */
double kron(int x, int y);

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

/**
 * @brief Computes the two-variable function H2(x, y).
 * @param x First input.
 * @param y Second input.
 * @return Result of H2.
 */
double H2(double x, double y);

/**
 * @brief Computes the Passarino-Veltman scalar function B.
 * @param m1 First mass.
 * @param m2 Second mass.
 * @param Q Energy scale.
 * @return Result of the B function.
 */
double B(double m1, double m2, double Q);

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

//for wilson coefficients

/**
 * @brief Scalar one-loop Wilson coefficient A0t.
 * @param x Input parameter.
 * @return Coefficient A0t(x).
 */
double A0t(double x);

/**
 * @brief Scalar one-loop Wilson coefficient F0t.
 * @param x Input parameter.
 * @return Coefficient F0t(x).
 */
double F0t(double x);

/**
 * @brief Scalar one-loop Wilson coefficient B0t.
 * @param x Input parameter.
 * @return Coefficient B0t(x).
 */
double B0t(double x);

/**
 * @brief Scalar one-loop Wilson coefficient C0t.
 * @param x Input parameter.
 * @return Coefficient C0t(x).
 */
double C0t(double x);

/**
 * @brief Scalar one-loop Wilson coefficient D0t.
 * @param x Input parameter.
 * @return Coefficient D0t(x).
 */
double D0t(double x);

/**
 * @brief Scalar one-loop Wilson coefficient E0t.
 * @param x Input parameter.
 * @return Coefficient E0t(x).
 */
double E0t(double x);

/**
 * @brief Wilson coefficient T(x).
 * @param x Input parameter.
 * @return Coefficient T(x).
 */
double T(double x);

/**
 * @brief Wilson coefficient A1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient A1t(x, l).
 */
double A1t(double x, double l);

/**
 * @brief Wilson coefficient B1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient B1t(x, l).
 */
double B1t(double x, double l);

/**
 * @brief Wilson coefficient C1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient C1t(x, l).
 */
double C1t(double x, double l);

/**
 * @brief Wilson coefficient D1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient D1t(x, l).
 */
double D1t(double x, double l);

/**
 * @brief Wilson coefficient E1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient E1t(x, l).
 */
double E1t(double x, double l);

/**
 * @brief Wilson coefficient F1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient F1t(x, l).
 */
double F1t(double x,double l);

/**
 * @brief Wilson coefficient G1t depending on x and scale l.
 * @param x Input parameter.
 * @param l Logarithmic scale.
 * @return Coefficient G1t(x, l).
 */
double G1t(double x, double l);

/**
 * @brief Wilson coefficient C7t at 2m_t scale.
 * @param x Input parameter.
 * @return Coefficient C7t(2m_t).
 */
double C7t2mt(double x);

/**
 * @brief Wilson coefficient C7c at M_W scale.
 * @param x Input parameter.
 * @return Coefficient C7c(M_W).
 */
double C7c2MW(double x);

/**
 * @brief Wilson coefficient C8t at 2m_t scale.
 * @param x Input parameter.
 * @return Coefficient C8t(2m_t).
 */
double C8t2mt(double x);

/**
 * @brief Wilson coefficient C8c at M_W scale.
 * @param x Input parameter.
 * @return Coefficient C8c(M_W).
 */
double C8c2MW(double x);

/**
 * @brief Wilson coefficient F7_1.
 * @param x Input parameter.
 * @return Coefficient F7_1(x).
 */
double F7_1(double x);

/**
 * @brief Wilson coefficient F7_2.
 * @param x Input parameter.
 * @return Coefficient F7_2(x).
 */
double F7_2(double x);

/**
 * @brief Wilson coefficient F8_1.
 * @param x Input parameter.
 * @return Coefficient F8_1(x).
 */
double F8_1(double x);

/**
 * @brief Wilson coefficient F8_2.
 * @param x Input parameter.
 * @return Coefficient F8_2(x).
 */
double F8_2(double x);

/**
 * @brief Wilson coefficient G3H depending on x and lu.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient G3H(x, lu).
 */
double G3H(double x, double lu);

/**
 * @brief Wilson coefficient G4H depending on x and lu.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient G4H(x, lu).
 */
double G4H(double x, double lu);

/**
 * @brief Wilson coefficient G7H depending on x, lu, and ld.
 * @param x Input parameter.
 * @param lu Upper log scale.
 * @param ld Lower log scale.
 * @return Coefficient G7H(x, lu, ld).
 */
double G7H(double x, double lu, double ld);

/**
 * @brief Wilson coefficient G8H depending on x, lu, and ld.
 * @param x Input parameter.
 * @param lu Upper log scale.
 * @param ld Lower log scale.
 * @return Coefficient G8H(x, lu, ld).
 */
double G8H(double x, double lu, double ld);

/**
 * @brief Wilson coefficient EH depending on x and lu.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient EH(x, lu).
 */
double EH(double x, double lu);

/**
 * @brief Wilson coefficient D9H0 depending on x and lu.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient D9H0(x, lu).
 */
double D9H0(double x, double lu);

/**
 * @brief Wilson coefficient D9H1 depending on x, lu, and L.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @param L Additional scale or parameter.
 * @return Coefficient D9H1(x, lu, L).
 */
double D9H1(double x, double lu, double L);

/**
 * @brief Wilson coefficient Delta3H depending on x and lu.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient Delta3H(x, lu).
 */
double Delta3H(double x, double lu);

/**
 * @brief Wilson coefficient Delta4H depending on x and lu.
 * @param x Input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient Delta4H(x, lu).
 */
double Delta4H(double x, double lu);

/**
 * @brief Wilson coefficient Delta7H depending on x, lu, and ld.
 * @param x Input parameter.
 * @param lu Upper log scale.
 * @param ld Lower log scale.
 * @return Coefficient Delta7H(x, lu, ld).
 */
double Delta7H(double x, double lu, double ld);

/**
 * @brief Wilson coefficient Delta8H depending on x, lu, and ld.
 * @param x Input parameter.
 * @param lu Upper log scale.
 * @param ld Lower log scale.
 * @return Coefficient Delta8H(x, lu, ld).
 */
double Delta8H(double x, double lu, double ld);

/**
 * @brief Wilson coefficient C9llH0 depending on x, y, and lu.
 * @param x First input parameter.
 * @param y Second input parameter.
 * @param lu Logarithmic scale.
 * @return Coefficient C9llH0(x, y, lu).
 */
double C9llH0(double x, double y, double lu);

/**
 * @brief Wilson coefficient C9llH1 depending on x, y, lu, and L.
 * @param x First input parameter.
 * @param y Second input parameter.
 * @param lu Logarithmic scale.
 * @param L Additional scale or parameter.
 * @return Coefficient C9llH1(x, y, lu, L).
 */
double C9llH1(double x, double y, double lu, double L);

/**
 * @brief Wilson coefficient C10Wt at 2m_t scale.
 * @param x Input parameter.
 * @return Coefficient C10Wt(2m_t).
 */
double C10Wt2mt(double x);

/**
 * @brief Wilson coefficient C10Wc at M_W scale.
 * @param x Input parameter.
 * @return Coefficient C10Wc(M_W).
 */
double C10Wc2MW(double x);

/**
 * @brief Wilson coefficient C10Zt at 2m_t scale.
 * @param x Input parameter.
 * @return Coefficient C10Zt(2m_t).
 */
double C10Zt2mt(double x);

/**
 * @brief Wilson coefficient C10Z from triangle diagrams.
 * @param x Input parameter.
 * @return Coefficient C10Z (triangle contribution).
 */
double C10Z2tri(double x);

/**
 * @brief Wilson coefficient F0SP depending on xt.
 * @param xt Input parameter.
 * @return Coefficient F0SP(xt).
 */
double F0SP(double xt);

/**
 * @brief Wilson special funtion S0 depending on xt.
 * @param xt Input parameter.
 * @return Value of S0(xt).
 */
double S0(double x);

/**
 * @brief Wilson special funtion D0 depending on 4 parameters.
 * @param w Input parameter.
 * @param x Input parameter.
 * @param y Input parameter.
 * @param z Input parameter.
 * @return Value of D0(w,x,y,z).
 */
double D0(double w, double x, double y, double z);


/**
 * @brief Wilson special funtion D2p depending on 4 parameters.
 * @param w Input parameter.
 * @param x Input parameter.
 * @param y Input parameter.
 * @param z Input parameter.
 * @return Value of D2p(w,x,y,z).
 */
double D2p(double w, double x, double y, double z);


/**
 * @brief Wilson special funtion getDelta.
 * @param delta Input parameter.
 * @param Z Input parameter.
 * @param M Input parameter.
 * @param m_av Input parameter.
 * @param delta_LL Input parameter.
 * @param delta_LR Input parameter.
 * @param delta_RL Input parameter.
 * @param delta_RR Input parameter.
 * @return Value of getDelta.
 */
void getDelta(scalar_t delta[6][6],scalar_t Z[6][6],double M[6],double m_av,scalar_t delta_LL[3][3],scalar_t delta_LR[3][3],scalar_t delta_RL[3][3],scalar_t delta_RR[3][3]);

/**
 * @brief Wilson special funtion h3 depending on x.
 * @param x Input parameter.
 * @return Value of h3(x).
 */
double h3(double x);

/**
 * @brief Wilson special funtion h1 depending on x.
 * @param x Input parameter.
 * @return Value of h1(x).
 */
double h1(double x);

/**
 * @brief Wilson special funtion h4 depending on x and y.
 * @param x Input parameter.
 * @param y Input parameter.
 * @return Value of h4(x,y).
 */
double h4(double x, double y);

/**
 * @brief Wilson special funtion f depending on x.
 * @param x Input parameter.
 * @return Value of f(x).
 */
double f(double x);

double Y0(double xt);
double Y1(double xt, double mu, double mass_W);
double X0(double xt);
double X1(double xt, double mu, double mass_W);

/*
*   Numerical constants for flavor observables 
*/ 

constexpr double GAMMA = 0.5772156649015328;

// Digamma function for integer inputs
double psi(int n);


/*****************************************
 *                                       *
 *      Observable-related functions     *
 *                                       *
 *****************************************/

// TODO : Docstring

namespace BV {
    complex_t A_Seidel  (double s_hat, double L_b);
    complex_t B_Seidel  (double s_hat, double L_b);
    complex_t C_Seidel  (double s, double mu_b);

    complex_t f_17      (double s_hat, double L_b, double z, size_t max_pow=20);
    complex_t f_27      (double s_hat, double L_b, double z, size_t max_pow=20);
    complex_t f_27_u    (double s_hat, double L_b);
    complex_t f_19_PS   (double s_hat, double L_b, double z, size_t max_pow=20);
    complex_t f_19_1S   (double s_hat, double L_b, double z, size_t max_pow=20);
    complex_t f_19_u    (double s_hat, double L_b);
    complex_t f_29_PS   (double s_hat, double L_b, double z, size_t max_pow=20);
    complex_t f_29_1S   (double s_hat, double L_b, double z, size_t max_pow=20);
    complex_t f_29_u    (double s_hat, double L_b);
    complex_t f_87      (double s_hat, double L_b);
    complex_t f_89      (double s_hat);

    complex_t h(double s, double m_q, double mu_b);
    complex_t B_0(double s, double m_q);
    complex_t L_1(complex_t x);
    complex_t I_1(double u, double s_hat, double m_q_hat);
    complex_t G(double x_bar, double z);
    complex_t hard_kernel(double u, double z);
    complex_t G2(double z, double L_b); 
    complex_t G8(double L_b);

}; // namespace BV

namespace KP {
    complex_t F(double z);
    complex_t G(double z);
    complex_t H(double z, double r_P);
}; // namespace KP

#endif