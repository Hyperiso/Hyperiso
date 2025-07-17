#include <cmath>
#include <vector>
#include <functional>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "Math.h"

/*
    Integration helper functions
*/

double unwrap_lambda(double x, void *p) {
    auto fun = static_cast<std::function<double(double)>*>(p);
    return (*fun)(x);
}

double integrate(Integrand f, double l, double u, double prec) {
    double res, err;
    size_t max_intervals = 1000;
    gsl_function F;
    F.function = &unwrap_lambda;
    F.params = &f;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_intervals);
    gsl_integration_qag(&F, l, u, 0, prec, max_intervals, GSL_INTEG_GAUSS21, w, &res, &err);
    gsl_integration_workspace_free(w);

    return res;
}

double cauchy_principal_value(Integrand f, double l, double u, double pole, double prec) {
    double res, err;
    size_t max_intervals = 1000;
    gsl_function F;
    F.function = &unwrap_lambda;
    F.params = &f;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_intervals);
    gsl_integration_qawc(&F, l, u, pole, 0, prec, max_intervals, w, &res, &err);
    gsl_integration_workspace_free(w);

    return res;
}

double cauchy_principal_value_s(Integrand f, double l, double u, double pole, std::vector<double> s, double prec) {
    double I = 0;

    for (size_t k = 0; k < s.size() - 1; k++) {
        if (s[k] < pole && pole < s[k + 1]) {
            I += cauchy_principal_value(f, s[k], s[k + 1], pole, prec);
        } else {
            I += integrate(f, s[k], s[k + 1], prec);
        }
    }

    return I;
}

scalar_t c_integrate(cIntegrand f, double l, double u, double prec) {
    return scalar_t(integrate([f] (double u) -> double { return f(u).real(); }, l, u, prec),
                     integrate([f] (double u) -> double { return f(u).imag(); }, l, u, prec));
}
