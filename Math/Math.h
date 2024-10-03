#pragma once
#include <cmath>
#include <complex>
#include "Formfactors.h"
#include <map>
#include <vector>

typedef std::complex<double> complex_t;

/* Constants */

constexpr double PI =       3.1415926535897932;
constexpr double PI2 =      9.869604401089358;
constexpr double PI3 =      31.006276680299816;
constexpr double INV_PI =   0.318309886183791;
constexpr double INV_PI2 =  0.101321183642338;
constexpr double INV_PI3 =  0.032251534433199;
constexpr double E =        2.7182818284590452;
constexpr double ZETA3 =    1.2020569031595942;
constexpr double RT2 =      1.414213562373095;
constexpr double INV_RT2 =  0.707106781186547;

constexpr double HBAR =     6.58211889e-25;  // GeV.s
constexpr double G_NEWT =   6.67428e-8;  // cm^3.g^-1.s^-2
constexpr double M_PLANCK = 1.2209102930946623e+19; // GeV
constexpr double M_NUCL =   0.939;  // Gev
constexpr double M_P =   0.9382720;  // Gev
constexpr double M_N =   0.9395654;  // Gev

constexpr std::complex<double> I = std::complex<double>(0, 1);

/* Functions */

double kron(int x, int y);

double Li2(double x);
complex_t CLi2(complex_t x);
double Cl2(double x);
double H2(double x, double y);
double B(double m1, double m2, double Q);
template <typename T> int sgn(T val) {
    return (T(0) < val) -(val<T(0));
}

//for wilson coefficients
double A0t(double x);
double F0t(double x);
double B0t(double x);
double C0t(double x);
double D0t(double x);
double E0t(double x);
double T(double x);

double A1t(double x, double l);
double B1t(double x, double l);
double C1t(double x, double l);
double D1t(double x, double l);
double E1t(double x, double l);
double F1t(double x,double l);
double G1t(double x, double l);

double C7t2mt(double x);
double C7c2MW(double x);
double C8t2mt(double x);
double C8c2MW(double x);

double F7_1(double x);
double F7_2(double x);
double F8_1(double x);
double F8_2(double x);

double G3H(double x, double lu);
double G4H(double x, double lu);
double G7H(double x, double lu, double ld);
double G8H(double x, double lu, double ld);

double EH(double x, double lu);


double D9H0(double x, double lu);
double D9H1(double x, double lu, double L);

double Delta3H(double x, double lu);
double Delta4H(double x, double lu);
double Delta7H(double x, double lu, double ld);
double Delta8H(double x, double lu, double ld);


double C9llH0(double x, double y, double lu);
double C9llH1(double x, double y, double lu, double L);

double C10Wt2mt(double x);
double C10Wc2MW(double x);
double C10Zt2mt(double x);
double C10Z2tri(double x);

double F0SP(double xt);

/*
    For observables
*/ 

// B > Ks gamma isospin asymmetry

complex_t G(double, double);
double F_perp(double a_1_perp, double a_2_perp); // Done
double X_perp(double a_1_perp, double a_2_perp, double cutoff); // Done
complex_t G_perp(double s, double a_1_perp, double a_2_perp); // Done
complex_t G2_perp(double s, double rb); // Done
complex_t G8_perp(double rb); // Done
complex_t H_perp(double s, double a_1_par, double a_2_par); // Done
complex_t H2_perp(double s, double a_1_perp, double a_2_perp); // Done
double H8_perp(double a_1_perp, double a_2_perp); // Done

using SparseMatrix = std::map<std::pair<std::string, std::string>, double>;

std::vector<std::string> getDiagonalElements(const SparseMatrix& matrix);
SparseMatrix invertMatrix(const SparseMatrix& matrix, const std::vector<std::string>& indices);
