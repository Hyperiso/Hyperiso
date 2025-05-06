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

scalar_t c_integrate(cIntegrand f, double l, double u, double prec) {
    return scalar_t(integrate([f] (double u) { return f(u).real(); }, l, u, prec),
                     integrate([f] (double u) { return f(u).imag(); }, l, u, prec));
}
