#include <cmath>
#include <complex>

#include <complex>

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
    if (abs(x) > 1) return -(hpl1(0, 1.0 / x) * (pi * i * -1 * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0)) - 2.0 * (hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + pi * i * 0.5 * pow(hpl1(0, 1.0 / x), 2) + pow(hpl1(0, 1.0 / x), 3) / 6.0);
    if (real(x) > 0.5) return -(hpl1(1, 1.0 - x) * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + pow(pi, 2) / 6.0)) - 2.0 * (1.2020569031595942 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - (hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2)) / 2.0);
    return hpl_base1(0, x) * hpl_base2(0, 1, x) - 2.0 * hpl_base3(0, 0, 1, x);
}

  if (i1 == 0 && i2 == 1 && i3 == 1) {
    if (abs(x) > 1) return 1.2020569031595942 + pi * i * -1 * hpl2(0, 1, 1.0 / x) - hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + i * 0.16666666666666666 * pow(pi, 3) + pi * i * -0.5 * pow(hpl1(0, 1.0 / x), 2) - pow(hpl1(0, 1.0 / x), 3) / 6.0;
    if (real(x) > 0.5) return 1.2020569031595942 + hpl1(0, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 0, 1, 1.0 - x) - (hpl1(1, 1.0 - x) * pow(hpl1(0, 1.0 - x), 2)) / 2.0;
    return hpl_base3(0, 1, 1, x);
}
  if (i1 == 1 && i2 == 0 && i3 == 0) {
    if (abs(x) > 1) return hpl3(0, 0, 1, 1.0 / x) - (hpl1(0, 1.0 / x) * pow(pi, 2)) / 3.0 + hpl1(0, 1.0 / x) * (pi * i * -1 * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) + pi * i * 0.5 * pow(hpl1(0, 1.0 / x), 2) + ((pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * pow(hpl1(0, 1.0 / x), 2)) / 2.0 + pow(hpl1(0, 1.0 / x), 3) / 6.0;
    if (real(x) > 0.5) return 1.2020569031595942 + hpl1(1, 1.0 - x) * hpl2(0, 1, 1.0 - x) - hpl3(0, 1, 1, 1.0 - x) + hpl1(1, 1.0 - x) * (hpl1(0, 1.0 - x) * hpl1(1, 1.0 - x) - hpl2(0, 1, 1.0 - x) + pow(pi, 2) / 6.0) - (hpl1(1, 1.0 - x) * pow(pi, 2)) / 6.0 - hpl1(0, 1.0 - x) * pow(hpl1(1, 1.0 - x), 2);
    return -(hpl_base1(0, x) * hpl_base2(0, 1, x)) + hpl_base3(0, 0, 1, x) + (hpl_base1(1, x) * pow(hpl_base1(0, x), 2)) / 2.0;
}

  if (i1 == 1 && i2 == 0 && i3 == 1) {
    if (abs(x) > 1) return (pi * i + hpl1(0, 1.0 / x) + hpl1(1, 1.0 / x)) * (pi * i * -1 * hpl1(0, 1.0 / x) - hpl2(0, 1, 1.0 / x) + pow(pi, 2) / 3.0 - pow(hpl1(0, 1.0 / x), 2) / 2.0) - 2.0 * (1.2020569031595942 + pi * i * -1 * hpl2(0, 1, 1.0 / x) - hpl1(0, 1.0 / x) * hpl2(0, 1, 1.0 / x) + hpl3(0, 0, 1, 1.0 / x) - hpl3(0, 1, 1, 1.0 / x) + (hpl1(0, 1.0 / x) * pow(pi, 2)) / 2.0 + i * 0.16666666666666666 * pow(pi, 3) + pi * i * -0.5 * pow(hpl1(0, 1.0 / x), 2) - pow(hpl1(0, 1.0 / x), 3) / 6.0);
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
