#include "integration.h"

#include <gsl/gsl_errno.h>

namespace {
    std::once_flag gsl_error_handler_once_flag;

    void ensure_gsl_nonfatal() {
        std::call_once(gsl_error_handler_once_flag, []() {
            gsl_set_error_handler_off();
        });
    }
}

double integrate(RealValuedFunction f, double l, double u, double prec) {
    ensure_gsl_nonfatal();
    if (!std::isfinite(l) || !std::isfinite(u) || !std::isfinite(prec) || prec <= 0.0) {
        throw std::invalid_argument("integrate: invalid bounds or precision");
    }
    if (l == u) {
        return 0.0;
    }
    double res = 0.0;
    double err = 0.0;
    size_t max_intervals = 1000;
    gsl_function F;
    F.function = &unwrap_lambda_unidim;
    F.params = &f;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_intervals);
    if (!w) {
        throw std::runtime_error("integrate: gsl_integration_workspace_alloc failed");
    }
    const int status = gsl_integration_qag(
        &F, l, u,
        0.0, prec,
        max_intervals,
        GSL_INTEG_GAUSS21,
        w,
        &res, &err
    );

    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS || !std::isfinite(res)) {
        std::ostringstream oss;
        oss << "integrate: gsl_integration_qag failed"
            << " status=" << status
            << " (" << gsl_strerror(status) << ")"
            << ", l=" << l
            << ", u=" << u
            << ", prec=" << prec
            << ", res=" << res
            << ", err=" << err;
        throw std::runtime_error(oss.str());
    }

    return res;

    return res;
}

scalar_t c_integrate(ComplexValuedFunction f, double l, double u, double prec) {
    return scalar_t(integrate([f] (double x) -> double { return f(x).real(); }, l, u, prec),
                     integrate([f] (double x) -> double { return f(x).imag(); }, l, u, prec));
}