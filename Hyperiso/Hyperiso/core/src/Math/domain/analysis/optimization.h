#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>
#include "functions.h"
#include "gsl_wrappers.h"

struct MinimizationContext {
    std::size_t simplex_max_iter {500};
    std::size_t bfgs_max_iter {500};
    double switch_tol {1e-3};
    double final_tol {1e-5};
    double simplex_initial_step_size {0.1};
    double bfgs_line_search_tol {0.1};
    double bfgs_initial_step_size {0.01};
};

struct MinimizationResult {
    int status;
    std::vector<double> argmin;
    double min;
};

struct ScalarMinimizationContext {
    std::array<double, 2> bracket;
    double start;
    double tol {1e-3};
    std::size_t max_iter {1000};
    bool sorted {false};
};

struct ScalarMinimizationResult {
    int status;
    double argmin;
    double min;
};

// TODO : use GSL for consistency and make return type MinimizationResult
bool find_bracket(const RealValuedFunction& f, double x_min, double x_max, double &a, double &b, int n_samples = 100);   
double brent_root(const RealValuedFunction& f, double a, double b, double xtol = 1e-6, double ftol = 1e-6, int max_it = 100);

template<typename Container, typename U>
std::size_t bisect_array(const Container& X, const U& x_0) {
    auto it = std::lower_bound(X.begin(), X.end(), x_0);

    if (it == X.begin()) return 0;
    if (it == X.end())   return X.size() - 1;

    return it - X.begin(); // faster than std::distance
}

ScalarMinimizationResult minimize_scalar(RealValuedFunction f, const ScalarMinimizationContext& context);

MinimizationResult minimize_NM(RealValuedForm f, const std::vector<double>& x0, const std::vector<double>& scales, const MinimizationContext& context);
MinimizationResult minimize_BFGS(RealValuedForm f, const std::vector<double>& x0, const std::vector<double>& scales, const MinimizationContext& context);
MinimizationResult minimize_NM(ScaledForm f, const std::vector<double>& t0, const MinimizationContext& context);
MinimizationResult minimize_BFGS(ScaledForm f, const std::vector<double>& t0, const MinimizationContext& context);
MinimizationResult minimize_combined(RealValuedForm f, const std::vector<double>& x0, const std::vector<double>& scales, const MinimizationContext& context);

#endif // OPTIMIZATION_H
