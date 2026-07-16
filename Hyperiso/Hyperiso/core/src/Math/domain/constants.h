#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <complex>

// Mathematical
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
constexpr double GAMMA =        0.5772156649015328;
constexpr double ERF_INV_RT2 =  0.682689492137085;

constexpr std::complex<double> I = std::complex<double>(0, 1);

// Physical
constexpr double HBAR =     6.58211889e-25;  // GeV.s
constexpr double G_NEWT =   6.67428e-8;  // cm^3.g^-1.s^-2
constexpr double M_PLANCK = 1.2209102930946623e+19; // GeV
constexpr double M_NUCL =   0.939;  // Gev
constexpr double M_P =   0.9382720;  // Gev
constexpr double M_N =   0.9395654;  // Gev
constexpr double GEV_TO_INV_PS = 1.519268e12;
constexpr double GEV_TO_INV_S = 1.519268e24;

#endif // __CONSTANTS_H__
