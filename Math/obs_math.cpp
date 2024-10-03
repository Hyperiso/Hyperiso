#include <cmath>
#include <vector>
#include <functional>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "Math.h"

/*
    Integration helper functions
*/

using Integrand = std::function<double(double)>;
using cIntegrand = std::function<complex_t(double)>;

double unwrap_lambda(double x, void *p) {
    auto fun = static_cast<Integrand*>(p);
    return (*fun)(x);
}

double integrate(Integrand f, double l, double u) {
    double res, err;
    size_t n_eval = 1000;
    gsl_function F;
    F.function = &unwrap_lambda;
    F.params = &f;

    gsl_integration_qng(&F, l, u, 0, 1e-4, &res, &err, &n_eval);

    return res;
}

complex_t c_integrate(cIntegrand f, double l, double u) {
    return complex_t(integrate([f] (double u) { return std::real(f(u)); }, l, u),
                     integrate([f] (double u) { return std::imag(f(u)); }, l, u));
}

/*
*   RGE related functions
*/


/*
    B -> K* gamma isospin asymmetry
*/

// Auxiliary functions

complex_t h(double u, double s) {
    complex_t rt = std::sqrt((u - 4 * s + I * 1e-10) / u);
    return (4 * s * (CLi2(2. / (1. - rt)) + CLi2(2. / (1. + rt))) / u - 2.) / u;
}

complex_t g_2(double s) {
    double ls = std::log(s);
    double ls2 = ls * ls;
    double ls3 = ls * ls2;
    complex_t a_0 = -833 / 162. - 20 * PI * I / 27.;
    complex_t a_1 = 48 - 5 * PI2 - 36 * ZETA3 + I * (30 * PI - 2 * PI3) + (36 - 9 * PI2 + 6 * PI * I) * ls + (3. + 6 * PI * I) * ls2 + ls3;
    complex_t a_2 = 18 + 2 * PI2 - 2 * PI3 * I + (12 - 6 * PI2) * ls + 6 * PI * I * ls2 + ls3;
    complex_t a_3 = -9 - 14 * PI2 + 112 * PI * I + (182. - 48 * PI * I) * ls - 126 * ls2;
    return a_0 + 2. * (s * a_1 + s * s * a_2) / 9. + s * s * s * a_3 / 27. + 8 * PI2 * std::pow(s, 1.5) / 9;
}

double phi_perp(double u, double a_1_perp, double a_2_perp) {
    double ubar = 1 - u;
    double xi = u - ubar;
    return 6 * u * ubar * (1 + 3 * a_1_perp * xi + 1.5 * a_2_perp * (5 * xi * xi - 1));
}

double ga_perp(double u, double a_1_par, double a_2_par) {
    double ubar = 1 - u;
    double xi = u - ubar;
    double uubar = u * ubar;
    double ul = u * std::log(u);
    double ulbar = ubar * std::log(ubar);
    double a2 = 0.25 * a_2_par + 5 * FFInput::zeta_3_A * (1 - 0.1875 * FFInput::omega_10_A) / 3 + 8.75 * FFInput::zeta_3_V;

    return 6 * (uubar * (1 + a_1_par * xi + a2 * (5 * xi * xi - 1)) 
                + FFInput::delta_tilde_p * (3 * uubar + ul + ulbar) 
                + FFInput::delta_tilde_m * (ulbar - ul));
}

double gv_perp(double u, double a_1_par, double a_2_par) {
    double ubar = 1 - u;
    double xi = u - ubar;
    double xi2 = xi * xi;
    double xi3 = xi2 * xi;
    double xi4 = xi2 * xi2;
    double l = std::log(u);
    double lbar = std::log(ubar);
    double a = 0.75;
    double b = 1.5 * a_1_par;
    double c = 3 * a_2_par / 7 + 5 * FFInput::zeta_3_A;
    double d = 9 * a_2_par / 112 + 6.5625 * FFInput::zeta_3_V - 0.234375 * FFInput::zeta_3_A * FFInput::omega_10_A; 
    
    return a * (1 + xi2) + b * xi3 + c * (3 * xi2 - 1) + d * (3 - 30 * xi2 + 35 * xi4)
            + 1.5 * (FFInput::delta_tilde_p * (2 + l + lbar) + FFInput::delta_tilde_m * (2 * xi + lbar - l));
}

complex_t G(double s, double xbar) {
    double eps = 1e-10;

    auto iG = [s, xbar, eps] (double u) {
        double v = u * (1 - u);
        return -4 * v * std::log(s - v * xbar - I * eps);
    };

    return c_integrate(iG, 0, 1);
}

// Header functions

double F_perp(double a_1_perp, double a_2_perp) {
    auto iF_perp = [a_1_perp, a_2_perp] (double u) {
        return phi_perp(u, a_1_perp, a_2_perp) / (3 * (1 - u));
    };

    return integrate(iF_perp, 0, 1);
}

complex_t G_perp(double s, double a_1_perp, double a_2_perp) {
    auto iG_perp = [s, a_1_perp, a_2_perp] (double x) {
        double xbar = 1 - x;
        return phi_perp(x, a_1_perp, a_2_perp) * G(s, xbar) / (3 * xbar);
    };

    return c_integrate(iG_perp, 0, 1);
}

complex_t G2_perp(double s, double rb) {
    return -104 * std::log(rb) / 27 + g_2(s);
}

complex_t G8_perp(double rb) {
    complex_t g_8 = 11 / 3 - 2 * PI2 / 9 + 2. * I * PI / 3.;
    return 8 * std::log(rb) / 3 + g_8;
}

complex_t H_perp(double s, double a_1_par, double a_2_par) {
    auto iH_perp = [s, a_1_par, a_2_par] (double x) {
        double xbar = 1 - x;
        double dx = 1e-5;
        double gv = gv_perp(x, a_1_par, a_2_par);
        double dga_dx = (ga_perp(x + dx, a_1_par, a_2_par) - ga_perp(x, a_1_par, a_2_par)) / dx;
        return (gv - dga_dx / 4) * G(s, xbar);
    }; 

    return c_integrate(iH_perp, 0, 1);
}

complex_t H2_perp(double s, double a_1_perp, double a_2_perp) {
    auto iH2_perp = [s, a_1_perp, a_2_perp] (double x) {
        double xbar = 1 - x;
        return -h(xbar, s) * phi_perp(x, a_1_perp, a_2_perp);
    }; 

    return c_integrate(iH2_perp, 0, 1);
}

double H8_perp(double a_1_perp, double a_2_perp) {
    auto iH8_perp = [a_1_perp, a_2_perp] (double x) {
        return phi_perp(x, a_1_perp, a_2_perp) / x;
    };

    return integrate(iH8_perp, 0, 1);
}

double X_perp(double a_1_perp, double a_2_perp, double cutoff) {
    auto iX_perp = [a_1_perp, a_2_perp] (double x) {
        double xbar = 1 - x;
        return phi_perp(x, a_1_perp, a_2_perp) * (1 + xbar) / (3 * xbar * xbar);
    };

    return integrate(iX_perp, 0, cutoff);
}