#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include "Math.h"

std::complex<double> cd(double x, double y) {
    return {x, y};
}

std::complex<double> hpl_base1(int i, std::complex<double> x) {
    if (i == 0) return std::log(x);
    if (i == 1) return -std::log(1.0 - x);
    return {0.0, 0.0};
}

std::complex<double> hpl_base2(int i1, int i2, std::complex<double> x) {
    std::complex<double> u = std::log(1.0 - x);

    if (i1 == 0 && i2 == 1) {

        return -u - 0.25 * std::pow(u, 2) - 0.027777777777777776 * std::pow(u, 3) + 
               0.0002777777777777778 * std::pow(u, 5) - 
               4.72411186696901e-6 * std::pow(u, 7) + 
               9.185773074661964e-8 * std::pow(u, 9) - 
               1.8978869988971e-9 * std::pow(u, 11) + 
               4.0647616451442256e-11 * std::pow(u, 13) - 
               8.921691020456452e-13 * std::pow(u, 15) + 
               1.9939295860721074e-14 * std::pow(u, 17) - 
               4.518980029619918e-16 * std::pow(u, 19) + 
               1.0356517612181247e-17 * std::pow(u, 21);
    }
    
    return {0.0, 0.0};
}


std::complex<double> hpl_base3(int i1, int i2, int i3, std::complex<double> x) {
    std::complex<double> u = std::log(1.0 - x);

    if (i1 == 0 && i2 == 0 && i3 == 1) {

        return -u - 0.375 * std::pow(u, 2) - 0.0787037037037037 * std::pow(u, 3) - 
               0.008680555555555556 * std::pow(u, 4) - 
               0.00012962962962962963 * std::pow(u, 5) + 
               0.00008101851851851852 * std::pow(u, 6) + 
               3.4193571608537595e-6 * std::pow(u, 7) - 
               1.328656462585034e-6 * std::pow(u, 8) - 
               8.660871756109851e-8 * std::pow(u, 9) + 
               2.52608759553204e-8 * std::pow(u, 10) + 
               2.144694468364065e-9 * std::pow(u, 11) - 
               5.140110622012979e-10 * std::pow(u, 12) - 
               5.24958211460083e-11 * std::pow(u, 13) + 
               1.0887754406636318e-11 * std::pow(u, 14) + 
               1.2779396094493695e-12 * std::pow(u, 15) - 
               2.369824177308745e-13 * std::pow(u, 16) - 
               3.104357887965462e-14 * std::pow(u, 17) + 
               5.261758629912506e-15 * std::pow(u, 18) + 
               7.538479549949265e-16 * std::pow(u, 19) - 
               1.1862322577752286e-16 * std::pow(u, 20) - 
               1.8316979965491384e-17 * std::pow(u, 21);
    }

    if (i1 == 0 && i2 == 1 && i3 == 1) {
        return 0.25 * std::pow(u, 2) + 0.08333333333333333 * std::pow(u, 3) + 
               0.010416666666666666 * std::pow(u, 4) - 
               0.00011574074074074075 * std::pow(u, 6) + 
               2.066798941798942e-6 * std::pow(u, 8) - 
               4.1335978835978836e-8 * std::pow(u, 10) + 
               8.698648744945042e-10 * std::pow(u, 12) - 
               1.887210763816962e-11 * std::pow(u, 14) + 
               4.182042665838962e-13 * std::pow(u, 16) - 
               9.415778600896063e-15 * std::pow(u, 18) + 
               2.146515514069461e-16 * std::pow(u, 20);
    }

    return {0.0, 0.0};
}

std::complex<double> hpl_base4(int i1, int i2, int i3, int i4, std::complex<double> x) {
    std::complex<double> u = std::log(1.0 - x);

    if (i1 == 0 && i2 == 0 && i3 == 0 && i4 == 1) {
        return -1.0 * u - 0.4375 * std::pow(u, 2) - 0.11651234567901235 * std::pow(u, 3) - 
               0.019820601851851853 * std::pow(u, 4) - 0.001927932098765432 * std::pow(u, 5) - 
               0.000031057098765432096 * std::pow(u, 6) + 0.000015624009114857836 * std::pow(u, 7) + 
               8.485123546773206e-7 * std::pow(u, 8) - 2.290961660318971e-7 * std::pow(u, 9) - 
               2.1832614218526917e-8 * std::pow(u, 10) + 3.882824879172015e-9 * std::pow(u, 11) + 
               5.446292103220332e-10 * std::pow(u, 12) - 6.960805210682725e-11 * std::pow(u, 13) - 
               1.3375737686445216e-11 * std::pow(u, 14) + 1.2784852685266572e-12 * std::pow(u, 15) + 
               3.260562858024892e-13 * std::pow(u, 16) - 2.364757116861826e-14 * std::pow(u, 17) - 
               7.923135122031162e-15 * std::pow(u, 18) + 4.3452915709984186e-16 * std::pow(u, 19) + 
               1.923627006253592e-16 * std::pow(u, 20) - 7.812414333195955e-18 * std::pow(u, 21);
    }

    if (i1 == 0 && i2 == 1 && i3 == 0 && i4 == 1) {
        return 0.25 * std::pow(u, 2) + 0.1111111111111111 * std::pow(u, 3) + 
               0.022569444444444444 * std::pow(u, 4) + 0.0020833333333333333 * std::pow(u, 5) - 
               0.000027006172839506174 * std::pow(u, 6) - 0.00001984126984126984 * std::pow(u, 7) + 
               4.527273872511968e-7 * std::pow(u, 8) + 3.389987682315725e-7 * std::pow(u, 9) - 
               7.939132443100697e-9 * std::pow(u, 10) - 6.6805622361177916e-9 * std::pow(u, 11) + 
               1.4490216610627064e-10 * std::pow(u, 12) + 1.39908336457158e-10 * std::pow(u, 13) - 
               2.7425719106565973e-12 * std::pow(u, 14) - 3.032441227329819e-12 * std::pow(u, 15) + 
               5.358569182999823e-14 * std::pow(u, 16) + 6.724068599976371e-14 * std::pow(u, 17) - 
               1.0756816626218996e-15 * std::pow(u, 18) - 1.5158529016922455e-15 * std::pow(u, 19) + 
               2.208955024323606e-17 * std::pow(u, 20) + 3.460964863954937e-17 * std::pow(u, 21);
    }

    if (i1 == 0 && i2 == 1 && i3 == 1 && i4 == 1) {
        return -0.05555555555555555 * std::pow(u, 3) - 0.020833333333333332 * std::pow(u, 4) - 
               0.002777777777777778 * std::pow(u, 5) + 0.00003306878306878307 * std::pow(u, 7) - 
               6.123848716441309e-7 * std::pow(u, 9) + 1.252605419272086e-8 * std::pow(u, 11) - 
               2.6765073061369356e-10 * std::pow(u, 13) + 5.871322376319437e-12 * std::pow(u, 15) - 
               1.312013385361243e-13 * std::pow(u, 17) + 2.97340376870402e-15 * std::pow(u, 19) - 
               6.814334965299877e-17 * std::pow(u, 21);
    }

    return {0.0, 0.0};
}

std::complex<double> hpl1(int i, std::complex<double> x) {
    if (i == 0) return std::log(x);
    if (i == 1) return -std::log(1.0 - x);
    return {0.0, 0.0};
}

std::complex<double> hpl2(int i1, int i2, std::complex<double> x) {
    const std::complex<double> pi = std::acos(-1);
    const std::complex<double> i(0.0, 1.0);

    if (i1 == 0 && i2 == 0) {
        if (std::abs(x) > 1) return std::pow(hpl1(0, 1.0 / x), 2) / 2.0;
        if (std::real(x) > 0.5) return std::pow(hpl1(1, 1.0 - x), 2) / 2.0;
        return std::pow(hpl_base1(0, x), 2) / 2.0;
    }

    if (i1 == 0 && i2 == 1) {
        if (std::abs(x) > 1) return
            pi * i * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + std::pow(pi, 2) / 3.0 - 
            std::pow(hpl1(0, 1.0 / x), 2) / 2.0;

        if (std::real(x) > 0.5) return
            hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + std::pow(pi, 2) / 6.0;

        return hpl_base2(0, 1, x);
    }

    if (i1 == 1 && i2 == 0) {
        if (std::abs(x) > 1) return
            pi * -i * hpl1(0, 1.0 / x) - 
            hpl1(0, 1.0 / x) * (pi * -i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) + 
            hpl2(0, 1, 1.0 / x) - std::pow(pi, 2) / 3.0 + std::pow(hpl1(0, 1.0 / x), 2) / 2.0;

        if (std::real(x) > 0.5) return hpl2(0, 1, 1.0 - x) - std::pow(pi, 2) / 6.0;

        return hpl_base1(0, x) * hpl_base1(1, x) - hpl_base2(0, 1, x);
    }

    if (i1 == 1 && i2 == 1) {
        if (std::abs(x) > 1) return std::pow(pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x), 2) / 2.0;
        if (std::real(x) > 0.5) return std::pow(hpl1(0, 1.0 - x), 2) / 2.0;
        return std::pow(hpl_base1(1, x), 2) / 2.0;
    }

    return {0.0, 0.0};
}

std::complex<double> hpl3(int i1, int i2, int i3, std::complex<double> x) {
    const std::complex<double> pi = std::acos(-1);
    const std::complex<double> i(0.0, 1.0);
    const double zeta3 = 1.2020569031595942; // Apery's constant

    if (i1 == 0 && i2 == 0 && i3 == 0) {
        if (std::abs(x) > 1) return -std::pow(hpl1(0, 1.0 / x), 3) / 6.0;
        if (std::real(x) > 0.5) return -std::pow(hpl1(1, 1.0 - x), 3) / 6.0;
        return std::pow(hpl_base1(0, x), 3) / 6.0;
    }

    if (i1 == 0 && i2 == 0 && i3 == 1) {
        if (std::abs(x) > 1) return hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * std::pow(pi, 2)) / 3.0 + pi * i * 0.5 * std::pow(hpl1(0, 1.0 / x), 2) + std::pow(hpl1(0, 1.0 / x), 3) / 6.0;
        if (std::real(x) > 0.5) return zeta3 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * std::pow(pi, 2)) / 6.0 - (hpl1(0, 1.0 - x) * std::pow(hpl1(1, 1.0 - x), 2)) / 2.0;
        return hpl_base3(0, 0, 1, x);
    }

  if (i1 == 0 && i2 == 1 && i3 == 0) {
    if (abs(x) > 1) return -(hpl1(0, 1.0 / x) * (-pi * i * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0)) - 2.0 * (hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + pi * i * 0.5 * pow(hpl1(0, 1.0 / x), 2) + pow(hpl1(0, 1.0 / x), 3) / 6.0);
    if (real(x) > 0.5) return -(hpl1(1, 1.0 - x) * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + pow(pi, 2) / 6.0)) - 2.0 * (1.2020569031595942 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0);
    return hpl_base1(0, x) * hpl_base2(0, 1, x) - 2.0 * hpl_base3(0, 0, 1, x);
}

  if (i1 == 0 && i2 == 1 && i3 == 1) {
    if (abs(x) > 1) return 1.2020569031595942 - pi * i * hpl2(0, 1, 1.0 / x) - hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + i * 0.16666666666666666 * pow(pi, 3) + pi * i * -0.5 * pow(hpl1(0, 1.0 / x), 2) - pow(hpl1(0, 1.0 / x), 3) / 6.0;
    if (real(x) > 0.5) return 1.2020569031595942 + hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0;
    return hpl_base3(0, 1, 1, x);
}
  if (i1 == 1 && i2 == 0 && i3 == 0) {
    if (abs(x) > 1) return hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + hpl1(0, 1.0 / x) * (-pi * i * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) + pi * i * 0.5 * pow(hpl1(0, 1.0 / x), 2) + ((pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + pow(hpl1(0, 1.0 / x), 3) / 6.0;
    if (real(x) > 0.5) return 1.2020569031595942 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) + hpl1(1, 1.0 - x) * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + pow(pi, 2) / 6.0) - (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2);
    return -(hpl_base1(0, x) * hpl_base2(0, 1, x)) + hpl_base3(0, 0, 1, x) + (hpl_base1(1, x) * pow(hpl_base1(0, x), 2)) / 2.0;
}

  if (i1 == 1 && i2 == 0 && i3 == 1) {
    if (abs(x) > 1) return (pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * (pi * i * -1. * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) - 2.0 * (1.2020569031595942 + pi * i * -1. * hpl2(0, 1, 1.0 / x) - hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + i * 0.16666666666666666 * pow(pi, 3) + pi * i * -0.5 * pow(hpl1(0, 1.0 / x), 2) - pow(hpl1(0, 1.0 / x), 3) / 6.0);
    if (real(x) > 0.5) return -(hpl1(0, 1.0 - x) * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + pow(pi, 2) / 6.0)) - 2.0 * (1.2020569031595942 + hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0);
    return hpl_base3(0, 0, 1, x) - (hpl_base1(0, x) * pow(pi, 2)) / 3.0 + pi * i * 0.5 * pow(hpl_base1(0, x), 2) + pow(hpl_base1(0, x), 3) / 6.0;
}
  if (i1 == 1 && i2 == 1 && i3 == 0) {
    if (abs(x) > 1) return -hpl3(0, 1, 1, 1.0 / x) + hpl2(0, 1, 1.0 / x) * (pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) - hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 - pi * i * 0.5 * pow(hpl1(0, 1.0 / x), 2) - pow(hpl1(0, 1.0 / x), 3) / 6.0;
    if (real(x) > 0.5) return -hpl3(0, 1, 1, 1.0 - x) + hpl2(0, 1, 1.0 - x) * hpl1(1, 1.0 - x) - hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) + (hpl1(0, 1.0 - x) * pow(pi, 2)) / 6.0 - (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0;
    return hpl_base3(0, 1, 0, x);
}

  if (i1 == 1 && i2 == 1 && i3 == 1) {
    if (abs(x) > 1) return -hpl3(0, 1, 1, 1.0 / x) + (pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * hpl2(0, 1, 1.0 / x) - hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + i * 0.16666666666666666 * pow(pi, 3) + pi * i * -0.5 * pow(hpl1(0, 1.0 / x), 2) - pow(hpl1(0, 1.0 / x), 3) / 6.0;
    if (real(x) > 0.5) return -hpl3(0, 1, 1, 1.0 - x) + hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0 + 1.2020569031595942;
    return hpl_base3(0, 1, 1, x);
}

    return {0.0, 0.0}; // Default return for unhandled cases
}

std::complex<double> hpl4(int i1, int i2, int i3, int i4, std::complex<double> x) {
    const std::complex<double> i(0, 1);
    const double pi = acos(-1);

    if (i1 == 0 && i2 == 0 && i3 == 0 && i4 == 0) {
        if (abs(x) > 1) return pow(hpl1(0, 1.0 / x), 4) / 24.0;
        if (real(x) > 0.5) return pow(hpl1(1, 1.0 - x), 4) / 24.0;
        return pow(hpl_base1(0, x), 4) / 24.0;
    }

    if (i1 == 0 && i2 == 0 && i3 == 0 && i4 == 1) {
        if (abs(x) > 1) return -hpl4(0, 0, 0, 1, 1.0 / x) + pow(pi, 4) / 45.0 + (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + pi * i * -0.16666666666666666 * pow(hpl1(0, 1.0 / x), 3) - pow(hpl1(0, 1.0 / x), 4) / 24.0;
        if (real(x) > 0.5) return -1.2020569031595942 * hpl1(1, 1.0 - x) + hpl1(1, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) - hpl4(0, 1, 1, 1, 1.0 - x) + pow(pi, 4) / 90.0 - (hpl2(0, 1, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0 + (pow(pi, 2) * pow(hpl1(1, 1.0 - x), 2)) / 12.0 + (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 3)) / 6.0;
        return hpl_base4(0, 0, 0, 1, x);
    }

    if (i1 == 0 && i2 == 0 && i3 == 1 && i4 == 0) {
        if (abs(x) > 1) return -(hpl1(0, 1.0 / x) * (hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + pi * i * 0.5 * pow(hpl1(0, 1.0 / x), 2) + pow(hpl1(0, 1.0 / x), 3) / 6.0)) - 3. * (-hpl4(0, 0, 0, 1, 1.0 / x) + pow(pi, 4) / 45.0 + (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + pi * i * -0.16666666666666666 * pow(hpl1(0, 1.0 / x), 3) - pow(hpl1(0, 1.0 / x), 4) / 24.0);
        if (real(x) > 0.5) return -(hpl1(1, 1.0 - x) * (1.2020569031595942 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0)) - 3. * (-1.2020569031595942 * hpl1(1, 1.0 - x) + hpl1(1, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) - hpl4(0, 1, 1, 1, 1.0 - x) + pow(pi, 4) / 90.0 - (hpl2(0, 1, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0 + (pow(pi, 2) * pow(hpl1(1, 1.0 - x), 2)) / 12.0 + (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 3)) / 6.0);
        return hpl_base1(0, x) * hpl_base3(0, 0, 1, x) - 3. * hpl_base4(0, 0, 0, 1, x);
    }

    if (i1 == 0 && i2 == 0 && i3 == 1 && i4 == 1) {
        if (abs(x) > 1) {
            return (-2. * (3.7763731361630786 * std::complex<double>(0, 2) + 
                          2.4041138063191885 * hpl1(0, 1.0 / x) + 
                          pi * std::complex<double>(0, 1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                          pi * std::complex<double>(0, -2) * hpl3(0, 0, 1, 1.0 / x) - 
                          2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                          4. * hpl4(0, 0, 0, 1, 1.0 / x) + hpl4(0, 1, 0, 1, 1.0 / x) - 
                          (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 + pow(pi, 4) / 90.0 + 
                          (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 - 
                          (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                          pi * std::complex<double>(0, 0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) + 
                          pow(hpl1(0, 1.0 / x), 4) / 24.0) + 
                   pow(pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                       pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0, 2)) / 4.0;
        } else if (real(x) > 0.5) {
            return (-2. * (2.4041138063191885 * hpl1(1, 1.0 - x) + 
               hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
               hpl4(0, 1, 0, 1, 1.0 - x) + 
               (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
               (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 + 
               pow(pi, 4) / 120.0 + 
               pow(hpl2(0, 1, 1.0 - x), 2) - 
               2. * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
               (2. * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
               (-2. * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) - 
               2. * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
               (2. * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
               (-2. * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0)) + 
               pow(hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
               pow(pi, 2) / 6.0, 2)) / 4.0;
        } else {
            return (-2. * hpl_base4(0, 1, 0, 1, x) + pow(hpl_base2(0, 1, x), 2)) / 4.0;
        }
    }

    if (i1 == 0 && i2 == 1 && i3 == 0 && i4 == 0) {
    if (abs(x) > 1) {
        return (((pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - 
                  pow(hpl1(0, 1.0 / x), 2) / 2.0) * pow(hpl1(0, 1.0 / x), 2) + 
                 4.0 * hpl1(0, 1.0 / x) * (hpl3(0, 0, 1, 1.0 / x) - 
                 (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + 
                 pi * std::complex<double>(0, 0.5) * pow(hpl1(0, 1.0 / x), 2) + 
                 pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                 6.0 * (-hpl4(0, 0, 0, 1, 1.0 / x) + pow(pi, 4) / 45.0 + 
                 (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                 pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                 pow(hpl1(0, 1.0 / x), 4) / 24.0)) / 2.0);
    } else if (real(x) > 0.5) {
        return ((hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + pow(pi, 2) / 6.0) * pow(hpl1(1, 1.0 - x), 2) + 
                4.0 * hpl1(1, 1.0 - x) * (1.2020569031595942 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - 
                (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0) + 
                6.0 * (-1.2020569031595942 * hpl1(1, 1.0 - x) + hpl1(1, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) - 
                hpl4(0, 1, 1, 1, 1.0 - x) + pow(pi, 4) / 90.0 - (hpl2(0, 1, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0 + 
                (pow(pi, 2) * pow(hpl1(1, 1.0 - x), 2)) / 12.0 + (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 3)) / 6.0)) / 2.0;
    } else {
        
        return 
		(-4.*hpl_base1(0,x)*hpl_base3(0,0,1,x) + 6.*hpl_base4(0,0,0,1,x) + 
		hpl_base2(0,1,x)*pow(hpl_base1(0,x),2))/2.;
    }
}
    if (i1 == 0 && i2 == 1 && i3 == 0 && i4 == 1) {
    if (abs(x) > 1) {
        return 3.7763731361630786 * std::complex<double>(0, 2) + 
               2.4041138063191885 * hpl1(0, 1.0 / x) + 
               pi * std::complex<double>(0, 1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
               pi * std::complex<double>(0, -2) * hpl3(0, 0, 1, 1.0 / x) - 
               2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
               4.0 * hpl4(0, 0, 0, 1, 1.0 / x) + 
               hpl4(0, 1, 0, 1, 1.0 / x) - 
               (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 + 
               pow(pi, 4) / 90.0 + 
               (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 - 
               (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
               pi * std::complex<double>(0, 0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) + 
               pow(hpl1(0, 1.0 / x), 4) / 24.0;
    } else if (real(x) > 0.5) {
        return 2.4041138063191885 * hpl1(1, 1.0 - x) + 
               hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
               hpl4(0, 1, 0, 1, 1.0 - x) + 
               (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
               (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 + 
               pow(pi, 4) / 120.0 + 
               pow(hpl2(0, 1, 1.0 - x), 2) - 
               2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
               (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
               (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) - 
               2.0 * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
               (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
               (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0);
    } else {
        return hpl_base4(0,1,0,1,x);
}
    }

    if (i1 == 0 && i2 == 1 && i3 == 1 && i4 == 0) {
        if (abs(x) > 1) {
            return 3.7763731361630786 * std::complex<double>(0, -2) - 
                   2.4041138063191885 * hpl1(0, 1.0 / x) + 
                   pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, 2) * hpl3(0, 0, 1, 1.0 / x) + 
                   2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) - 
                   4. * hpl4(0, 0, 0, 1, 1.0 / x) - 
                   hpl4(0, 1, 0, 1, 1.0 / x) + 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 - 
                   pow(pi, 4) / 90.0 - 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 - 
                   hpl1(0, 1.0 / x) * (1.2020569031595942 + 
                   pi * std::complex<double>(0, -1) * hpl2(0, 1, 1.0 / x) - 
                   hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   hpl3(0, 0, 1, 1.0 / x) - 
                   hpl3(0, 1, 1, 1.0 / x) + 
                   (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                   std::complex<double>(0, 0.16666666666666666) * pow(pi, 3) + 
                   pi * std::complex<double>(0, -0.5) * pow(hpl1(0, 1.0 / x), 2) - 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                   pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                   (2.0 * (3.7763731361630786 * std::complex<double>(0, 2) + 
                   2.4041138063191885 * hpl1(0, 1.0 / x) + 
                   pi * std::complex<double>(0, 1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, -2) * hpl3(0, 0, 1, 1.0 / x) - 
                   2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                   4. * hpl4(0, 0, 0, 1, 1.0 / x) + 
                   hpl4(0, 1, 0, 1, 1.0 / x) - 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 + 
                   pow(pi, 4) / 90.0 + 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 - 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                   pi * std::complex<double>(0, 0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) + 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0) - 
                   pow(pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                   pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0, 2)) / 2.0;
        } else if (real(x) > 0.5) {
            return -2.4041138063191885 * hpl1(1, 1.0 - x) - 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) + 
                   hpl4(0, 1, 0, 1, 1.0 - x) - 
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 + 
                   (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   pow(pi, 4) / 120.0 - 
                   hpl1(1, 1.0 - x) * (1.2020569031595942 + 
                   hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
                   hpl3(0, 0, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0) - 
                   pow(hpl2(0, 1, 1.0 - x), 2) + 
                   2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) + 
                   2.0 * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) + 
                   (2.0 * (2.4041138063191885 * hpl1(1, 1.0 - x) + 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
                   hpl4(0, 1, 0, 1, 1.0 - x) + 
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 + 
                   pow(pi, 4) / 120.0 + 
                   pow(hpl2(0, 1, 1.0 - x), 2)) - 
                   2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) - 
                   2.0 * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) - 
                   pow(hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                   pow(pi, 2) / 6.0, 2)) / 2.0;
        }
        else {
            hpl_base1(0,x)*hpl_base3(0,1,1,x) - hpl_base4(0,1,0,1,x) + 
   		(2.*hpl_base4(0,1,0,1,x) - pow(hpl_base2(0,1,x),2))/2.;
        }
    }

    if (i1 == 0 && i2 == 1 && i3 == 1 && i4 == 1) {
        if (abs(x) > 1) {
            return pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, 1) * hpl3(0, 0, 1, 1.0 / x) + 
                   hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, -1) * hpl3(0, 1, 1, 1.0 / x) - 
                   hpl1(0, 1.0 / x) * hpl3(0, 1, 1, 1.0 / x) - 
                   hpl4(0, 0, 0, 1, 1.0 / x) - 
                   hpl4(0, 1, 0, 1, 1.0 / x) / 2.0 - 
                   hpl4(0, 1, 1, 1, 1.0 / x) + 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                   std::complex<double>(0, 0.16666666666666666) * hpl1(0, 1.0 / x) * pow(pi, 3) - 
                   (19 * pow(pi, 4)) / 360.0 - 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 4.0 + 
                   pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 / x) - pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0 + 
                   pow(hpl2(0, 1, 1.0 / x), 2) / 4.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 / x) + pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0;
        } else if (real(x) > 0.5) {
            return hpl1(0, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) - 
                   hpl4(0, 0, 0, 1, 1.0 - x) + 
                   pow(pi, 4) / 90.0 - 
                   (hpl2(0, 1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0 + 
                   (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 3)) / 6.0;
        } else {
            return hpl_base4(0, 1, 1, 1, x);
        }
    }

    if (i1 == 1 && i2 == 0 && i3 == 0 && i4 == 0) {
        if (abs(x) > 1) {
            return (-3. * (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                          pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) * pow(hpl1(0, 1.0 / x), 2) - 
                    6. * hpl1(0, 1.0 / x) * (hpl3(0, 0, 1, 1.0 / x) - 
                    (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + 
                    pi * std::complex<double>(0, 0.5) * pow(hpl1(0, 1.0 / x), 2) + 
                    pow(hpl1(0, 1.0 / x), 3) / 6.0) - 
                    (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * pow(hpl1(0, 1.0 / x), 3) - 
                    6. * (-hpl4(0, 0, 0, 1, 1.0 / x) + pow(pi, 4) / 45.0 + 
                    (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                    pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                    pow(hpl1(0, 1.0 / x), 4) / 24.0)) / 6.0;
        } else if (real(x) > 0.5) {
            return (-3. * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                          pow(pi, 2) / 6.0) * pow(hpl1(1, 1.0 - x), 2) - 
                    6. * hpl1(1, 1.0 - x) * (1.2020569031595942 + 
                    hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - 
                    (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                    (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0) + 
                    hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 3) - 
                    6. * (-1.2020569031595942 * hpl1(1, 1.0 - x) + 
                    hpl1(1, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) - 
                    hpl4(0, 1, 1, 1, 1.0 - x) + pow(pi, 4) / 90.0 - 
                    (hpl2(0, 1, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0 + 
                    (pow(pi, 2) * pow(hpl1(1, 1.0 - x), 2)) / 12.0 + 
                    (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 3)) / 6.0)) / 6.0;
        } else {
            return (6. * hpl_base1(0, x) * hpl_base3(0, 0, 1, x) - 6. * hpl_base4(0, 0, 0, 1, x) - 
                    3. * hpl_base2(0, 1, x) * pow(hpl_base1(0, x), 2) + hpl_base1(1, x) * pow(hpl_base1(0, x), 3)) / 6.0;
        }
    }

    if (i1 == 1 && i2 == 0 && i3 == 0 && i4 == 1) {
        if (abs(x) > 1) {
            return 3.7763731361630786 * std::complex<double>(0, -2) - 
                   2.4041138063191885 * hpl1(0, 1.0 / x) + 
                   pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, 2) * hpl3(0, 0, 1, 1.0 / x) + 
                   2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) - 
                   4. * hpl4(0, 0, 0, 1, 1.0 / x) - 
                   hpl4(0, 1, 0, 1, 1.0 / x) + 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 - 
                   pow(pi, 4) / 90.0 - 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                   (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * 
                   (hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + 
                   pi * std::complex<double>(0, 0.5) * pow(hpl1(0, 1.0 / x), 2) + 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                   pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                   (2.0 * (3.7763731361630786 * std::complex<double>(0, 2) + 
                   2.4041138063191885 * hpl1(0, 1.0 / x) + 
                   pi * std::complex<double>(0, 1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, -2) * hpl3(0, 0, 1, 1.0 / x) - 
                   2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                   4. * hpl4(0, 0, 0, 1, 1.0 / x) + hpl4(0, 1, 0, 1, 1.0 / x) - 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 + pow(pi, 4) / 90.0 + 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 - 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                   pi * std::complex<double>(0, 0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) + 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0) - 
                   pow(pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                   pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0, 2)) / 2.0;
        } else if (real(x) > 0.5) {
            return -2.4041138063191885 * hpl1(1, 1.0 - x) - 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) + 
                   hpl4(0, 1, 0, 1, 1.0 - x) - 
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 + 
                   (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 - pow(pi, 4) / 120.0 - 
                   hpl1(0, 1.0 - x) * (1.2020569031595942 + 
                   hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0) - 
                   pow(hpl2(0, 1, 1.0 - x), 2) + 
                   2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) + 
                   2.0 * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) + 
                   (2.0 * (2.4041138063191885 * hpl1(1, 1.0 - x) + 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
                   hpl4(0, 1, 0, 1, 1.0 - x) + 
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 + pow(pi, 4) / 120.0 + 
                   pow(hpl2(0, 1, 1.0 - x), 2)) - 
                   2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) - 
                   pow(hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                   pow(pi, 2) / 6.0, 2)) / 2.0;
        } else {
            return hpl_base1(1, x) * hpl_base3(0, 0, 1, x) - hpl_base4(0, 1, 0, 1, x) + 
                   (2.0 * hpl_base4(0, 1, 0, 1, x) - pow(hpl_base2(0, 1, x), 2)) / 2.0;
        }
    }

    if (i1 == 1 && i2 == 0 && i3 == 1 && i4 == 0) {
        if (abs(x) > 1) {
            return 3.7763731361630786 * std::complex<double>(0, -2) - 
                   2.4041138063191885 * hpl1(0, 1.0 / x) + 
                   pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, 2) * hpl3(0, 0, 1, 1.0 / x) + 
                   2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) - 
                   4. * hpl4(0, 0, 0, 1, 1.0 / x) - 
                   hpl4(0, 1, 0, 1, 1.0 / x) + 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 - 
                   pow(pi, 4) / 90.0 - hpl1(0, 1.0 / x) * 
                   (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * 
                   (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                   pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) - 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                   2.0 * hpl1(0, 1.0 / x) * (1.2020569031595942 + 
                   pi * std::complex<double>(0, -1) * hpl2(0, 1, 1.0 / x) - 
                   hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - 
                   hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                   std::complex<double>(0, 0.16666666666666666) * pow(pi, 3) + 
                   pi * std::complex<double>(0, -0.5) * pow(hpl1(0, 1.0 / x), 2) - 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) - 
                   2.0 * (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * 
                   (hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + 
                   pi * std::complex<double>(0, 0.5) * pow(hpl1(0, 1.0 / x), 2) + 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                   pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                   pow(pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                   pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0, 2);
        } else if (real(x) > 0.5) {
            return -2.4041138063191885 * hpl1(1, 1.0 - x) - 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) + 
                   hpl4(0, 1, 0, 1, 1.0 - x) + 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * 
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                   pow(pi, 2) / 6.0) - (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * 
                   pow(pi, 2)) / 6.0 + (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   pow(pi, 4) / 120.0 + 2. * hpl1(1, 1.0 - x) * 
                   (1.2020569031595942 + hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
                   hpl3(0, 0, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0) + 
                   2.0 * hpl1(0, 1.0 - x) * (1.2020569031595942 + 
                   hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0) - 
                   pow(hpl2(0, 1, 1.0 - x), 2) + 
                   2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) + 
                   2.0 * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) + 
                   pow(hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                   pow(pi, 2) / 6.0, 2);
        } else {
            return hpl_base1(0, x) * hpl_base1(1, x) * hpl_base2(0, 1, x) - 
                   2. * hpl_base1(1, x) * hpl_base3(0, 0, 1, x) - 
                   2.0 * hpl_base1(0, x) * hpl_base3(0, 1, 1, x) - 
                   hpl_base4(0, 1, 0, 1, x) + 
                   pow(hpl_base2(0, 1, x), 2);
        }
    }

    if (i1 == 1 && i2 == 0 && i3 == 1 && i4 == 1) {
        if (abs(x) > 1) {
            return (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) *
                   (1.2020569031595942 + pi * std::complex<double>(0, -1) * hpl2(0, 1, 1.0 / x) - 
                   hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - 
                   hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                   std::complex<double>(0, 0.16666666666666666) * pow(pi, 3) + 
                   pi * std::complex<double>(0, -0.5) * pow(hpl1(0, 1.0 / x), 2) - 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) - 
                   3. * (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, 1) * hpl3(0, 0, 1, 1.0 / x) + 
                   hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, -1) * hpl3(0, 1, 1, 1.0 / x) - 
                   hpl1(0, 1.0 / x) * hpl3(0, 1, 1, 1.0 / x) - hpl4(0, 0, 0, 1, 1.0 / x) - 
                   hpl4(0, 1, 0, 1, 1.0 / x) / 2.0 - hpl4(0, 1, 1, 1, 1.0 / x) + 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                   std::complex<double>(0, 0.16666666666666666) * hpl1(0, 1.0 / x) * pow(pi, 3) - 
                   (19 * pow(pi, 4)) / 360.0 - 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 4.0 + 
                   pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 / x) - pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0 + 
                   pow(hpl2(0, 1, 1.0 / x), 2) / 4.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 / x) + pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0);
        } else if (real(x) > 0.5) {
            return -(hpl1(0, 1.0 - x) * (1.2020569031595942 + 
                   hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0)) - 
                   3. * (hpl1(0, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) - 
                   hpl4(0, 0, 0, 1, 1.0 - x) + pow(pi, 4) / 90.0 - 
                   (hpl2(0, 1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0 + 
                   (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 3)) / 6.0);
        } else {
            return hpl_base1(1, x) * hpl_base3(0, 1, 1, x) - 3. * hpl_base4(0, 1, 1, 1, x);
        }
    }

    if (i1 == 1 && i2 == 1 && i3 == 0 && i4 == 0) {
        if (abs(x) > 1) {
            return hpl1(0, 1.0 / x) * (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) *
                   (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                   pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) - 
                   hpl1(0, 1.0 / x) * (1.2020569031595942 + 
                   pi * std::complex<double>(0, -1) * hpl2(0, 1, 1.0 / x) - 
                   hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - 
                   hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                   std::complex<double>(0, 0.16666666666666666) * pow(pi, 3) + 
                   pi * std::complex<double>(0, -0.5) * pow(hpl1(0, 1.0 / x), 2) - 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                   (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) *
                   (hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + 
                   pi * std::complex<double>(0, 0.5) * pow(hpl1(0, 1.0 / x), 2) + 
                   pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                   (pow(hpl1(0, 1.0 / x), 2) * pow(pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + 
                   hpl1(1, 1.0 / x), 2)) / 4.0 + 
                   (2.0 * (3.7763731361630786 * std::complex<double>(0, 2) + 
                   2.4041138063191885 * hpl1(0, 1.0 / x) + 
                   pi * std::complex<double>(0, 1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                   pi * std::complex<double>(0, -2) * hpl3(0, 0, 1, 1.0 / x) - 
                   2.0 * hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                   4. * hpl4(0, 0, 0, 1, 1.0 / x) + hpl4(0, 1, 0, 1, 1.0 / x) - 
                   (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 3.0 + pow(pi, 4) / 90.0 + 
                   (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 - 
                   (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 6.0 + 
                   pi * std::complex<double>(0, 0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) + 
                   pow(hpl1(0, 1.0 / x), 4) / 24.0) - 
                   pow(pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                   pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0, 2)) / 4.0;
        } else if (real(x) > 0.5) {
            return -(hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) *
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                   pow(pi, 2) / 6.0)) - 
                   hpl1(1, 1.0 - x) * (1.2020569031595942 + 
                   hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0) + 
                   (pow(hpl1(0, 1.0 - x), 2) * pow(hpl1(1, 1.0 - x), 2)) / 4.0 - 
                   hpl1(0, 1.0 - x) * (1.2020569031595942 + 
                   hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - 
                   (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0) + 
                   (2.0 * (2.4041138063191885 * hpl1(1, 1.0 - x) + 
                   hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - 
                   hpl4(0, 1, 0, 1, 1.0 - x) + 
                   (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - 
                   (hpl2(0, 1, 1.0 - x) * pow(pi, 2)) / 6.0 + pow(pi, 4) / 120.0 + 
                   pow(hpl2(0, 1, 1.0 - x), 2) - 
                   2.0 * (hpl1(1, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - 
                   pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + 
                   pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0) - 
                   2.0 * (hpl1(0, 1.0 - x) * hpl3(0, 1, 1, 1.0 - x) + 
                   (2.0 * hpl4(0, 1, 0, 1, 1.0 - x) - 
                   pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0 + 
                   (-2.0 * hpl4(0, 1, 0, 1, 1.0 - x) + 
                   pow(hpl2(0, 1, 1.0 - x), 2)) / 2.0)) - 
                   pow(hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                   pow(pi, 2) / 6.0, 2)) / 4.0;
        } else {
            return -(hpl_base1(0, x) * hpl_base1(1, x) * hpl_base2(0, 1, x)) + hpl_base1(1, x) * hpl_base3(0, 0, 1, x) +
                   hpl_base1(0, x) * hpl_base3(0, 1, 1, x) + 
                   (pow(hpl_base1(0, x), 2) * pow(hpl_base1(1, x), 2)) / 4.0 + 
                   (2.0 * hpl_base4(0, 1, 0, 1, x) - pow(hpl_base2(0, 1, x), 2)) / 4.0;
        }
    }

    if (i1 == 1 && i2 == 1 && i3 == 0 && i4 == 1) {
        if (abs(x) > 1) {
            return (-4. * (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) *
                    (1.2020569031595942 + pi * std::complex<double>(0, -1) * hpl2(0, 1, 1.0 / x) - 
                    hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - 
                    hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                    std::complex<double>(0, 0.16666666666666666) * pow(pi, 3) + 
                    pi * std::complex<double>(0, -0.5) * pow(hpl1(0, 1.0 / x), 2) - 
                    pow(hpl1(0, 1.0 / x), 3) / 6.0) + 
                    (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                    pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) *
                    pow(pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x), 2) + 
                    6. * (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                    pi * std::complex<double>(0, 1) * hpl3(0, 0, 1, 1.0 / x) + 
                    hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                    pi * std::complex<double>(0, -1) * hpl3(0, 1, 1, 1.0 / x) - 
                    hpl1(0, 1.0 / x) * hpl3(0, 1, 1, 1.0 / x) - hpl4(0, 0, 0, 1, 1.0 / x) - 
                    hpl4(0, 1, 0, 1, 1.0 / x) / 2.0 - hpl4(0, 1, 1, 1, 1.0 / x) + 
                    (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                    std::complex<double>(0, 0.16666666666666666) * hpl1(0, 1.0 / x) * pow(pi, 3) - 
                    (19 * pow(pi, 4)) / 360.0 - 
                    (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                    (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 4.0 + 
                    pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                    pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                    (2.0 * hpl4(0, 1, 0, 1, 1.0 / x) - pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0 + 
                    pow(hpl2(0, 1, 1.0 / x), 2) / 4.0 + 
                    (-2.0 * hpl4(0, 1, 0, 1, 1.0 / x) + pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0)) / 2.0;
        } else if (real(x) > 0.5) {
            return ((hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                    pow(pi, 2) / 6.0) * pow(hpl1(0, 1.0 - x), 2) + 
                    4. * hpl1(0, 1.0 - x) * (1.2020569031595942 + 
                    hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - 
                    (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0) + 
                    6. * (hpl1(0, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) - 
                    hpl4(0, 0, 0, 1, 1.0 - x) + pow(pi, 4) / 90.0 - 
                    (hpl2(0, 1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0 + 
                    (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 3)) / 6.0)) / 2.0;
        }
        return (-4. * hpl_base1(1, x) * hpl_base3(0, 1, 1, x) + 6. * hpl_base4(0, 1, 1, 1, x) + 
                hpl_base2(0, 1, x) * pow(hpl_base1(1, x), 2)) / 2.0;
    }


    if (i1 == 1 && i2 == 1 && i3 == 1 && i4 == 0) {
        if (abs(x) > 1) {
            return (6. * (pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) *
                    (1.2020569031595942 + pi * std::complex<double>(0, -1) * hpl2(0, 1, 1.0 / x) - 
                    hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - 
                    hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                    std::complex<double>(0, 0.16666666666666666) * pow(pi, 3) + 
                    pi * std::complex<double>(0, -0.5) * pow(hpl1(0, 1.0 / x), 2) - 
                    pow(hpl1(0, 1.0 / x), 3) / 6.0) - 
                    3. * (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + 
                    pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) *
                    pow(pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x), 2) - 
                    hpl1(0, 1.0 / x) * pow(pi * std::complex<double>(0, 1) + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x), 3) - 
                    6. * (pi * std::complex<double>(0, -1) * hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + 
                    pi * std::complex<double>(0, 1) * hpl3(0, 0, 1, 1.0 / x) + 
                    hpl1(0, 1.0 / x) * hpl3(0, 0, 1, 1.0 / x) + 
                    pi * std::complex<double>(0, -1) * hpl3(0, 1, 1, 1.0 / x) - 
                    hpl1(0, 1.0 / x) * hpl3(0, 1, 1, 1.0 / x) - hpl4(0, 0, 0, 1, 1.0 / x) - 
                    hpl4(0, 1, 0, 1, 1.0 / x) / 2.0 - hpl4(0, 1, 1, 1, 1.0 / x) + 
                    (hpl2(0, 1, 1.0 / x) * pow(pi, 2)) / 2.0 + 
                    std::complex<double>(0, 0.16666666666666666) * hpl1(0, 1.0 / x) * pow(pi, 3) - 
                    (19 * pow(pi, 4)) / 360.0 - 
                    (hpl2(0, 1, 1.0 / x) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + 
                    (pow(pi, 2) * pow(hpl1(0, 1.0 / x), 2)) / 4.0 + 
                    pi * std::complex<double>(0, -0.16666666666666666) * pow(hpl1(0, 1.0 / x), 3) - 
                    pow(hpl1(0, 1.0 / x), 4) / 24.0 + 
                    (2.0 * hpl4(0, 1, 0, 1, 1.0 / x) - pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0 + 
                    pow(hpl2(0, 1, 1.0 / x), 2) / 4.0 + 
                    (-2.0 * hpl4(0, 1, 0, 1, 1.0 / x) + pow(hpl2(0, 1, 1.0 / x), 2)) / 2.0)) / 6.0;
        } else if (real(x) > 0.5) {
            return (-3. * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + 
                    pow(pi, 2) / 6.0) * pow(hpl1(0, 1.0 - x), 2) - 
                    6. * hpl1(0, 1.0 - x) * (1.2020569031595942 + 
                    hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - 
                    (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0) + 
                    hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 3) - 
                    6. * (hpl1(0, 1.0 - x) * hpl3(0, 0, 1, 1.0 - x) - 
                    hpl4(0, 0, 0, 1, 1.0 - x) + pow(pi, 4) / 90.0 - 
                    (hpl2(0, 1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0 + 
                    (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 3)) / 6.0)) / 6.0;
        }
        return (6. * hpl_base1(1, x) * hpl_base3(0, 1, 1, x) - 6. * hpl_base4(0, 1, 1, 1, x) - 
                3. * hpl_base2(0, 1, x) * pow(hpl_base1(1, x), 2) + hpl_base1(0, x) * pow(hpl_base1(1, x), 3)) / 6.0;
    }


    if(i1==1&&i2==1&&i3==1&&i4==1)
	{
		if(abs(x)>1) return pow(pi*cd(0.,1.) + hpl1(0,1./x) + hpl1(1,1./x),4)/24.;
		if(real(x)>.5) return pow(hpl1(0,1. - x),4)/24.;
		return pow(hpl_base1(1,x),4)/24.;
	}  


	return {0.,0.};

}


template <size_t MAX_N, size_t MAX_M, size_t MAX_LE, size_t MAX_NC>
void initialize_constants(std::array<std::array<double, MAX_M>, MAX_N> &s1, std::array<std::array<double, MAX_M>, MAX_N> &c, std::array<std::array<double, MAX_LE>, MAX_NC> &a) {
    // Initialisation des constantes s1
    s1[1][1] = 1.6449340668482;
    s1[1][2] = 1.2020569031596;
    s1[1][3] = 1.0823232337111;
    s1[1][4] = 1.0369277551434;
    s1[2][1] = 1.2020569031596;
    s1[2][2] = 2.7058080842778e-1;
    s1[2][3] = 9.6551159989444e-2;
    s1[3][1] = 1.0823232337111;
    s1[3][2] = 9.6551159989444e-2;
    s1[4][1] = 1.0369277551434;

    // Initialisation des constantes c
    c[1][1] = 1.6449340668482;
    c[1][2] = 1.2020569031596;
    c[1][3] = 1.0823232337111;
    c[1][4] = 1.0369277551434;
    c[2][1] = 0.0;
    c[2][2] = -1.8940656589945;
    c[2][3] = -3.0142321054407;
    c[3][1] = 1.8940656589945;
    c[3][2] = 3.0142321054407;
    c[4][1] = 0.0;



    a[0][1]=.96753215043498;
    a[1][1]=.16607303292785;
    a[2][1]=.02487932292423;
    a[3][1]=.00468636195945;
    a[4][1]=.00100162749616;
    a[5][1]=.00023200219609;
    a[6][1]=.00005681782272;
    a[7][1]=.00001449630056;
    a[8][1]=.00000381632946;
    a[9][1]=.00000102990426;
    a[10][1]=.00000028357538;
    a[11][1]=.00000007938705;
    a[12][1]=.00000002253670;
    a[13][1]=.00000000647434;
    a[14][1]=.00000000187912;
    a[15][1]=.00000000055029;
    a[16][1]=.00000000016242;
    a[17][1]=.00000000004827;
    a[18][1]=.00000000001444;
    a[19][1]=.00000000000434;
    a[20][1]=.00000000000131;
    a[21][1]=.00000000000040;
    a[22][1]=.00000000000012;
    a[23][1]=.00000000000004;
    a[24][1]=.00000000000001;
    a[0][2]=.95180889127832;
    a[1][2]=.43131131846532;
    a[2][2]=.10002250714905;
    a[3][2]=.02442415595220;
    a[4][2]=.00622512463724;
    a[5][2]=.00164078831235;
    a[6][2]=.00044407920265;
    a[7][2]=.00012277494168;
    a[8][2]=.00003453981284;
    a[9][2]=.00000985869565;
    a[10][2]=.00000284856995;
    a[11][2]=.00000083170847;
    a[12][2]=.00000024503950;
    a[13][2]=.00000007276496;
    a[14][2]=.00000002175802;
    a[15][2]=.00000000654616;
    a[16][2]=.00000000198033;
    a[17][2]=.00000000060204;
    a[18][2]=.00000000018385;
    a[19][2]=.00000000005637;
    a[20][2]=.00000000001735;
    a[21][2]=.00000000000536;
    a[22][2]=.00000000000166;
    a[23][2]=.00000000000052;
    a[24][2]=.00000000000016;
    a[25][2]=.00000000000005;
    a[26][2]=.00000000000002;
    a[0][3]=.98161027991365;
    a[1][3]=.72926806320726;
    a[2][3]=.22774714909321;
    a[3][3]=.06809083296197;
    a[4][3]=.02013701183064;
    a[5][3]=.00595478480197;
    a[6][3]=.00176769013959;
    a[7][3]=.00052748218502;
    a[8][3]=.00015827461460;
    a[9][3]=.00004774922076;
    a[10][3]=.00001447920408;
    a[11][3]=.00000441154886;
    a[12][3]=.00000135003870;
    a[13][3]=.00000041481779;
    a[14][3]=.00000012793307;
    a[15][3]=.00000003959070;
    a[16][3]=.00000001229055;
    a[17][3]=.00000000382658;
    a[18][3]=.00000000119459;
    a[19][3]=.00000000037386;
    a[20][3]=.00000000011727;
    a[21][3]=.00000000003687;
    a[22][3]=.00000000001161;
    a[23][3]=.00000000000366;
    a[24][3]=.00000000000116;
    a[25][3]=.00000000000037;
    a[26][3]=.00000000000012;
    a[27][3]=.00000000000004;
    a[28][3]=.00000000000001;
    a[0][4]=1.0640521184614;
    a[1][4]=1.0691720744981;
    a[2][4]=.41527193251768;
    a[3][4]=.14610332936222;
    a[4][4]=.04904732648784;
    a[5][4]=.01606340860396;
    a[6][4]=.00518889350790;
    a[7][4]=.00166298717324;
    a[8][4]=.00053058279969;
    a[9][4]=.00016887029251;
    a[10][4]=.00005368328059;
    a[11][4]=.00001705923313;
    a[12][4]=.00000542174374;
    a[13][4]=.00000172394082;
    a[14][4]=.00000054853275;
    a[15][4]=.00000017467795;
    a[16][4]=.00000005567550;
    a[17][4]=.00000001776234;
    a[18][4]=.00000000567224;
    a[19][4]=.00000000181313;
    a[20][4]=.00000000058012;
    a[21][4]=.00000000018579;
    a[22][4]=.00000000005955;
    a[23][4]=.00000000001911;
    a[24][4]=.00000000000614;
    a[25][4]=.00000000000197;
    a[26][4]=.00000000000063;
    a[27][4]=.00000000000020;
    a[28][4]=.00000000000007;
    a[29][4]=.00000000000002;
    a[30][4]=.00000000000001;
    a[0][5]=.97920860669175;
    a[1][5]=.08518813148683;
    a[2][5]=.00855985222013;
    a[3][5]=.00121177214413;
    a[4][5]=.00020722768531;
    a[5][5]=.00003996958691;
    a[6][5]=.00000838064065;
    a[7][5]=.00000186848945;
    a[8][5]=.00000043666087;
    a[9][5]=.00000010591733;
    a[10][5]=.00000002647892;
    a[11][5]=.00000000678700;
    a[12][5]=.00000000177654;
    a[13][5]=.00000000047342;
    a[14][5]=.00000000012812;
    a[15][5]=.00000000003514;
    a[16][5]=.00000000000975;
    a[17][5]=.00000000000274;
    a[18][5]=.00000000000077;
    a[19][5]=.00000000000022;
    a[20][5]=.00000000000006;
    a[21][5]=.00000000000002;
    a[22][5]=.00000000000001;
    a[0][6]=.95021851963952;
    a[1][6]=.29052529161433;
    a[2][6]=.05081774061716;
    a[3][6]=.00995543767280;
    a[4][6]=.00211733895031;
    a[5][6]=.00047859470550;
    a[6][6]=.00011334321308;
    a[7][6]=.00002784733104;
    a[8][6]=.00000704788108;
    a[9][6]=.00000182788740;
    a[10][6]=.00000048387492;
    a[11][6]=.00000013033842;
    a[12][6]=.00000003563769;
    a[13][6]=.00000000987174;
    a[14][6]=.00000000276586;
    a[15][6]=.00000000078279;
    a[16][6]=.00000000022354;
    a[17][6]=.00000000006435;
    a[18][6]=.00000000001866;
    a[19][6]=.00000000000545;
    a[20][6]=.00000000000160;
    a[21][6]=.00000000000047;
    a[22][6]=.00000000000014;
    a[23][6]=.00000000000004;
    a[24][6]=.00000000000001;
    a[0][7]=.95064032186777;
    a[1][7]=.54138285465171;
    a[2][7]=.13649979590321;
    a[3][7]=.03417942328207;
    a[4][7]=.00869027883583;
    a[5][7]=.00225284084155;
    a[6][7]=.00059516089806;
    a[7][7]=.00015995617766;
    a[8][7]=.00004365213096;
    a[9][7]=.00001207474688;
    a[10][7]=.00000338018176;
    a[11][7]=.00000095632476;
    a[12][7]=.00000027313129;
    a[13][7]=.00000007866968;
    a[14][7]=.00000002283195;
    a[15][7]=.00000000667205;
    a[16][7]=.00000000196191;
    a[17][7]=.00000000058018;
    a[18][7]=.00000000017246;
    a[19][7]=.00000000005151;
    a[20][7]=.00000000001545;
    a[21][7]=.00000000000465;
    a[22][7]=.00000000000141;
    a[23][7]=.00000000000043;
    a[24][7]=.00000000000013;
    a[25][7]=.00000000000004;
    a[26][7]=.00000000000001;
    a[0][8]=.98800011672229;
    a[1][8]=.04364067609601;
    a[2][8]=.00295091178278;
    a[3][8]=.00031477809720;
    a[4][8]=.00004314846029;
    a[5][8]=.00000693818230;
    a[6][8]=.00000124640350;
    a[7][8]=.00000024293628;
    a[8][8]=.00000005040827;
    a[9][8]=.00000001099075;
    a[10][8]=.00000000249467;
    a[11][8]=.00000000058540;
    a[12][8]=.00000000014127;
    a[13][8]=.00000000003492;
    a[14][8]=.00000000000881;
    a[15][8]=.00000000000226;
    a[16][8]=.00000000000059;
    a[17][8]=.00000000000016;
    a[18][8]=.00000000000004;
    a[19][8]=.00000000000001;
    a[0][9]=.95768506546350;
    a[1][9]=.19725249679534;
    a[2][9]=.02603370313918;
    a[3][9]=.00409382168261;
    a[4][9]=.00072681707110;
    a[5][9]=.00014091879261;
    a[6][9]=.00002920458914;
    a[7][9]=.00000637631144;
    a[8][9]=.00000145167850;
    a[9][9]=.00000034205281;
    a[10][9]=.00000008294302;
    a[11][9]=.00000002060784;
    a[12][9]=.00000000522823;
    a[13][9]=.00000000135066;
    a[14][9]=.00000000035451;
    a[15][9]=.00000000009436;
    a[16][9]=.00000000002543;
    a[17][9]=.00000000000693;
    a[18][9]=.00000000000191;
    a[19][9]=.00000000000053;
    a[20][9]=.00000000000015;
    a[21][9]=.00000000000004;
    a[22][9]=.00000000000001;
    a[0][10]=.99343651671347;
    a[1][10]=.02225770126826;
    a[2][10]=.00101475574703;
    a[3][10]=.00008175156250;
    a[4][10]=.00000899973547;
    a[5][10]=.00000120823987;
    a[6][10]=.00000018616913;
    a[7][10]=.00000003174723;
    a[8][10]=.00000000585215;
    a[9][10]=.00000000114739;
    a[10][10]=.00000000023652;
    a[11][10]=.00000000005082;
    a[12][10]=.00000000001131;
    a[13][10]=.00000000000259;
    a[14][10]=.00000000000061;
    a[15][10]=.00000000000015;
    a[16][10]=.00000000000004;
    a[17][10]=.00000000000001;
}






std::complex<double> polylog(int n, int m, double x) {

    const int MAX_N = 5;
    const int MAX_M = 5;
    const int MAX_LE = 31;
    const int MAX_NC = 11;

    std::array<std::array<double, MAX_M>, MAX_N> s1, c;
    std::array<std::array<double, MAX_LE>, MAX_NC> a;

    initialize_constants<MAX_N, MAX_M, MAX_LE, MAX_NC>(s1, c, a);


    double u[5],s1[5][5],c[5][5],a[31][11];
    const int fct[5] = {1, 1, 2, 6, 24};
    const int sgn[5] = {1, -1, 1, -1, 1};
    const int index[32] = {0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 8, 9, 0, 0, 0, 0, 0, 0, 0, 0, 10};
    const int nc[11] = {0, 24, 26, 28, 30, 22, 24, 26, 19, 22, 17};

    std::complex<double> z = 0.0;
    std::complex<double> v[6] = {0.0};
    std::complex<double> sk;
    std::complex<double> sj;

    double z1=1.;
    double hf=z1/2.;
    double b0{0.};

    if ((n < 1) || (n > 4) || (m < 1) || (m > 4) || ((n + m) > 5)) {
        return 0.0;
    }

    if (x == 1.0) {
        z = s1[n][m];
    } else if ((x > 2.0) || (x < -1.0)) {
        // Logique pour x > 2 ou x < -1
        double x1 = 1.0 / x;
        double h = 4.0 / 3.0 * x1 + 1.0 / 3.0;
        double alfa = h + h;
        v[0] = 1.0;
        v[1] = log(-x);

        for (int le = 2; le <= n + m; le++) {
            v[le] = v[1] * v[le - 1] / (le*1.);
        }


        for (int ke = 0; ke <= m - 1; ke++) {
            int m1 = m - ke;
            double r = pow(x1, m1) / ((fct[m1] * fct[n - 1])*1.);
            sj = 0.0;

            for (int je = 0; je <= ke; je++) {
                int n1 = n + ke - je;
                int le = index[10 * n1 + m1 - 10];
                double b1 = 0.0;
                double b2 = 0.0;

                for (int it = nc[le]; it >= 0; it--) {
                    b0 = a[it][le] + alfa * b1 - b2;
                    b2 = b1;
                    b1 = b0;
                }

                double q = (fct[n1 - 1] / fct[ke - je]) * (b0 - h * b2) * r / pow(m1, n1);
                sj += v[je] * q;
            }

            sk += sgn[ke]*1. * sj;
        }



        for (int je = 0; je <= n - 1; je++) {
            sj += v[je] * c[n - je][m];
        }

        std::complex<double> z = sgn[n] *1.* sk + sgn[m]*1. * (sj + v[n + m]);
    } else if (x>hf) {
        // Logique pour x > 0.5
        double x1 = 1.0 - x;
        double h = 4.0 / 3.0 * x1 + 1.0 / 3.0;
        double alfa = h + h;
        v[0] = 1.0;
        u[0] = 1.0;
        v[1] = log(x1);
        u[1] = log(x);

        for (int le = 2; le <= m; le++) {
            v[le] = v[1] * v[le - 1] / (1.*le);
        }

        for (int le = 2; le <= n; le++) {
            u[le] = u[1] * u[le - 1] / le;
        }

        sk = 0.0;

        for (int ke = 0; ke <= n - 1; ke++) {
            int m1 = n - ke;
            double r = pow(x1, m1) / fct[m1];

            for (int je = 0; je <= m - 1; je++) {
                int n1 = m - je;
                int le = index[10 * n1 + m1 - 10];
                double b1 = 0.0;
                double b2 = 0.0;

                for (int it = nc[le]; it >= 0; it--) {
                    b0 = a[it][le] + alfa * b1 - b2;
                    b2 = b1;
                    b1 = b0;
                }

                double q = sgn[je] * (b0 - h * b2) * r / pow(m1, n1);
                sj += v[je] * q;
            }

            sk = sk + u[ke] * (s1[m1][m] - sj);
        }

        z = sk + sgn[m] * u[n] * v[m];
    } else {
        // Logique pour les autres cas
        int le = index[10 * n + m - 10];
        double h = 4.0 / 3.0 * x + 1.0 / 3.0;
        double alfa = h + h;
        double b1 = 0.0;
        double b2 = 0.0;

        for (int it = nc[le]; it >= 0; it--) {
            b0 = a[it][le] + alfa * b1 - b2;
            b2 = b1;
            b1 = b0;
        }

        z = (b0 - h * b2) * pow(x, m) / (fct[m] * pow(m, n));
    }

    return z;
}

// LI2 ! -------------------------------------------------------------------------------------------------------------------

double Li2(double x) {
    const double pisq6 = M_PI * M_PI / 6.0;
    const double x_0 = -0.3;
    const double x_1 = 0.25;
    const double x_2 = 0.51;

    if (x == 1.0) {
        return pisq6;
    }

    if (x <= x_0) {
        double temp = log(fabs(1.0 - x));
        return -Li2(-x / (1.0 - x)) - temp * temp / 2.0;
    } else if (x < x_1) {
        double z = -log(1.0 - x);
        double temp = z * (1.0 - z / 4.0 * (1.0 - z / 9.0 * (1.0 - z * z / 100.0 * (1.0 - 5.0 * z * z / 294.0 * (1.0 - 7.0 * z * z / 360.0 * (1.0 - 5.0 * z * z / 242.0 * (1.0 - 7601.0 * z * z / 354900.0 * (1.0 - 91.0 * z * z / 4146.0 * (1.0 - 3617.0 * z * z / 161840.0)))))))));
        return temp;
    } else if (x < x_2) {
        return -Li2(-x) + Li2(x * x) / 2.0;
    } else {
        return pisq6 - Li2(1.0 - x) - log(fabs(x)) * log(fabs(1.0 - x));
    }
}


// LI3 ! -------------------------------------------------------------------------------------------------------------------


const double zeta3 = 1.202056903159594; // Apery's constant

double Li3(double x) {
    const double pisq6 = 16.0 * M_PI * M_PI / 6.0;
    const double x_0 = -1.0;
    const double x_1 = -0.85;
    const double x_2 = 0.25;
    const double x_3 = 0.63;
    const double x_4 = 1.0;

    if (x == 1.0) return zeta3;
    if (x == -1.0) return -0.75 * zeta3;

    if (x <= x_0) {
        double lnx = log(-x);
        return Li3(1.0 / x) - pisq6 * lnx - lnx * lnx * lnx / 6.0;
    } else if (x < x_1) {
        return Li3(x * x) / 4.0 - Li3(-x);
    } else if (x < x_2) {
    double z = -log(1.0 - x);
    double temp = z * (1.0 - 3.0 * z / 8.0 * (1.0 - 17.0 * z / 81.0 * (1.0 - 15.0 * z / 136.0 *
        (1.0 - 28.0 * z / 1875.0 * (1.0 + 5.0 * z / 8.0 * (1.0 - 304.0 * z / 7203.0 *
        (1.0 + 945.0 * z / 2432.0 * (1.0 - 44.0 * z / 675.0 * (1.0 + 7.0 * z / 24.0 *
        (1.0 - 26104.0 * z / 307461.0 * (1.0 + 1925.0 * z / 8023.0 *
        (1.0 - 53598548.0 * z / 524808375.0 *
        (1.0 + 22232925.0 * z / 107197096.0)))))))))))));
    return temp;
} else if (x < x_3) {
        return Li3(x * x) / 4.0 - Li3(-x);
    } else if (x < x_4) {
        double ln1x = log(1.0 - x);
        return -Li3(1.0 - x) - Li3(-x / (1.0 - x)) + zeta3 + pisq6 * ln1x - log(x) * ln1x * ln1x / 2.0 + ln1x * ln1x * ln1x / 6.0;
    } else {
        double lnx = log(x);
        return Li3(1.0 / x) + 2.0 * pisq6 * lnx - lnx * lnx * lnx / 6.0;
    }
}

std::complex<double> Li4(double x)
/* calculates the quadrilogarithm function of x */
{
	return polylog(3,1,x);
}

/* calculates the dilogarithm function of x, extended to complex numbers */ 
/*-------------------------------------------------------------*/

std::complex<double> CLi2(std::complex<double> x) {
    const double pisq6 = std::pow((4.0 * std::atan(1.0)), 2.0) / 6.0;

    const double x_0 = -0.30;
    const double x_1 = 0.25;
    const double x_2 = 0.51;

    if (x == 1.0) {
        return pisq6;
    }

    if (std::real(x) >= x_2) {
        return pisq6 - CLi2(1.0 - x) - std::log(x) * std::log(1.0 - x);
    }

    if ((std::abs(std::imag(x)) > 1.0) || (std::norm(x) > 1.2)) {
        return -CLi2(1.0 / x) - 0.5 * std::log(-x) * std::log(-x) - pisq6;
    }

    if (std::real(x) <= x_0) {
        std::complex<double> zz = std::log(1.0 - x);
        return -CLi2(-x / (1.0 - x)) - zz * zz / 2.0;
    } else if (std::real(x) < x_1) {
        std::complex<double> z = -std::log(1.0 - x);
        std::complex<double> temp = z * (1.0 - z / 4.0 * (1.0 - z / 9.0 * (1.0 - z * z / 100.0 *
            (1.0 - 5.0 * z * z / 294.0 * (1.0 - 7.0 * z * z / 360.0 *
            (1.0 - 5.0 * z * z / 242.0 * (1.0 - 7601.0 * z * z / 354900.0 *
            (1.0 - 91.0 * z * z / 4146.0 * (1.0 - 3617.0 * z * z / 161840.0)))))))));
        return temp;
    } else {
        return -CLi2(-x) + CLi2(x * x) / 2.0;
    }
}


std::complex<double> CLi3(std::complex<double> x)
/* calculates the trilogarithm function of x, extended to complex numbers */
{
	return hpl3(0,0,1,x);
}

std::complex<double> CLi4(std::complex<double> x)
/* calculates the quadrilogarithm function of x, extended to complex numbers */
{
	return hpl4(0,0,0,1,x);
}


double Cl2(double x) {
    std::complex<double> z = std::cos(x) + std::complex<double>(0, 1) * std::sin(x);
    return std::imag(CLi2(z));
}

double Cl3(double x) {
    std::complex<double> z = std::cos(x) + std::complex<double>(0, 1) * std::sin(x);
    return std::real(CLi3(z));
}




