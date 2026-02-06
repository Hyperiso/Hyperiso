#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <vector>
#include <functional>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>
#include "scalar.h"
#include "Matrix.h"

using RealValuedFunction = std::function<double(double)>;
using RealValuedForm = std::function<double(std::vector<double>)>;
using ComplexValuedFunction = std::function<scalar_t(double)>;

inline double unwrap_lambda_unidim(double x, void *p) {
    auto fun = static_cast<std::function<double(double)>*>(p);
    return (*fun)(x);
}

inline double unwrap_lambda_multidim(const gsl_vector* x, void *p) {
    auto fun = static_cast<std::function<double(std::vector<double>)>*>(p);
    std::vector<double> args (x->size);
    for (size_t i = 0; i < x->size; i++) {
        args[i] = gsl_vector_get(x, i);
    }
    return (*fun)(args);
}

/*
    Integration routines
*/

/**
 * @brief Performs numerical integration of a real-valued function of a real variable.
 * @param f Function to integrate.
 * @param l Lower bound of integration.
 * @param u Upper bound of integration.
 * @param prec Desired precision.
 * @return Numerical integral result.
 */
double integrate(RealValuedFunction f, double l, double u, double prec);

/**
 * @brief Performs numerical integration of a complex-valued function of a real variable.
 * @param f Complex function to integrate.
 * @param l Lower bound of integration.
 * @param u Upper bound of integration.
 * @param prec Desired precision.
 * @return Complex-valued result of the integral.
 */
scalar_t c_integrate(ComplexValuedFunction f, double l, double u, double prec);

/*
    Root finding and minimization algorithms
*/

bool find_bracket(const RealValuedFunction& f,
                  double x_min, double x_max,
                  double &a, double &b,
                  int n_samples = 100);
                  
/*
 * brent_root: Brent's method for root finding.
 *
 * f       : function<double(double)> to find root of (must be continuous)
 * a, b    : initial bracket with f(a) and f(b) of opposite signs
 * xtol    : absolute/relative tolerance stopping criterion for x
 * ftol    : tolerance for |f(x)|
 * max_it  : max iterations
 *
 * returns approximate root x (throws on invalid bracket or failure)
 */
double brent_root(const RealValuedFunction& f,
                  double a, double b,
                  double xtol = 1e-6,
                  double ftol = 1e-6,
                  int max_it = 100);

struct MinimizationContext {
    std::size_t max_iter {500};
    double tol {1e-5};
    std::vector<double> step_sizes;
};

struct MinimizationResult {
    int status;
    std::vector<double> argmin;
    double min;
};

MinimizationResult minimize(
    RealValuedForm f,
    const std::vector<double> x_start,
    const MinimizationContext& context
);

MinimizationResult minimize_scaled(
    RealValuedForm f,
    const std::vector<double>& x0,
    const MinimizationContext& context
);

/*
    Differentiation
*/

std::vector<double> gradient(RealValuedForm f, const std::vector<double>& x);
RealMatrix hessian(RealValuedForm f, const std::vector<double>& x);

#endif // __ANALYSIS_H__
