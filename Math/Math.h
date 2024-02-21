#pragma once
#include <cmath>

constexpr double pi = 3.141592654;


double Li2(double x);
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
