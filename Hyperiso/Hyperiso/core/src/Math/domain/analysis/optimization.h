#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <gsl/gsl_multimin.h>
#include "functions.h"
#include "gsl_wrappers.h"

struct MinimizationContext {
    std::size_t max_iter {500};
    double tol {1e-5};
    double step_size {0.1};
    double line_search_tol {0.1};
};

struct MinimizationResult {
    int status;
    std::vector<double> argmin;
    double min;
};

// TODO : use GSL for consistency and make return type MinimizationResult
bool find_bracket(const RealValuedFunction& f, double x_min, double x_max, double &a, double &b, int n_samples = 100);   
double brent_root(const RealValuedFunction& f, double a, double b, double xtol = 1e-6, double ftol = 1e-6, int max_it = 100);

MinimizationResult minimize_NM(RealValuedForm f, const std::vector<double> x_start, const MinimizationContext& context);
MinimizationResult minimize_BFGS(RealValuedForm f, const std::vector<double>& x0, const MinimizationContext& context);

#endif // OPTIMIZATION_H
