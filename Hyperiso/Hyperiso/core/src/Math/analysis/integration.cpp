#include "analysis.h"

/*
    Integration helper functions
*/

double integrate(RealValuedFunction f, double l, double u, double prec) {
    double res, err;
    size_t max_intervals = 1000;
    gsl_function F;
    F.function = &unwrap_lambda_unidim;
    F.params = &f;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_intervals);
    gsl_integration_qag(&F, l, u, 0, prec, max_intervals, GSL_INTEG_GAUSS21, w, &res, &err);
    gsl_integration_workspace_free(w);

    return res;
}

scalar_t c_integrate(ComplexValuedFunction f, double l, double u, double prec) {
    return scalar_t(integrate([f] (double u) -> double { return f(u).real(); }, l, u, prec),
                     integrate([f] (double u) -> double { return f(u).imag(); }, l, u, prec));
}
