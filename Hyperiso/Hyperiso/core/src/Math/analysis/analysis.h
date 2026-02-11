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

struct ScaledForm {
    RealValuedForm f;
    std::vector<double> x0;
    std::vector<double> s;

    ScaledForm(const RealValuedForm& f, const std::vector<double>& x0) : f(f), x0(x0) {
        this->s = std::vector(x0);
        std::transform(this->s.begin(), this->s.end(), this->s.begin(), std::abs<double>);
        double max_abs_x = *std::max_element(this->s.begin(), this->s.end());
        double eps = std::numeric_limits<double>::epsilon();
        double sigma_min = std::sqrt(eps) * (1 + max_abs_x);
        std::transform(this->s.begin(), this->s.end(), this->s.begin(), [sigma_min] (double x) { return std::max(x, sigma_min); });
    }

    double operator()(const std::vector<double>& t) const {
        std::vector<double> x(t.size());
        for (size_t i = 0; i < t.size(); i++) {
            x[i] = this->x0[i] + this->s[i] * t[i];

            std::cout << "t_i = " << t[i] << ", s_i = " << s[i] << ", x0_i = " << x0[i] << ", x_i = " << x[i] << std::endl;
        }

        return this->f(x);
    }   
};

/*
    Differentiation
*/

std::vector<double> gradient(const RealValuedForm& f, const std::vector<double>& x);
std::vector<double> gradient(const ScaledForm& f, const std::vector<double>& x);
RealMatrix hessian(const RealValuedForm& f, const std::vector<double>& x);
RealMatrix hessian(const ScaledForm& f, const std::vector<double>& t);
RealMatrix inverse_hessian(const RealValuedForm& f, const std::vector<double>& x);

/*
    GSL Wrappers
*/

inline double unwrap_lambda_unidim(double x, void *p) {
    auto fun = static_cast<std::function<double(double)>*>(p);
    return (*fun)(x);
}

inline double unwrap_lambda_multidim(const gsl_vector *t, void *params) {
    auto fun = static_cast<ScaledForm*>(params);
    std::vector<double> args (t->size);
    for (size_t i = 0; i < t->size; i++) {
        args[i] = gsl_vector_get(t, i);
    }
    return (*fun)(args);
}

inline void unwrap_lambda_gradient(const gsl_vector *t, void *params, gsl_vector *g) {
    auto fun = static_cast<ScaledForm*>(params);
    std::vector<double> args (t->size);
    for (size_t i = 0; i < t->size; i++) {
        args[i] = gsl_vector_get(t, i);
    }
    auto grad = gradient(*fun, args);
    for (size_t i = 0; i < grad.size(); i++) {
        gsl_vector_set(g, i, grad[i]);
    }
}

inline void unwrap_lambda_func_and_gradient(const gsl_vector *t, void *params, double *f, gsl_vector *g) {
    auto fun = static_cast<ScaledForm*>(params);
    std::vector<double> args (t->size);
    for (size_t i = 0; i < t->size; i++) {
        args[i] = gsl_vector_get(t, i);
    }
    *f = (*fun)(args);
    auto grad = gradient(*fun, args);
    for (size_t i = 0; i < grad.size(); i++) {
        gsl_vector_set(g, i, grad[i]);
    }
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
    double step_size {0.1};
    double line_search_tol {0.1};
};

struct MinimizationResult {
    int status;
    std::vector<double> argmin;
    double min;
};

MinimizationResult minimize_NM(
    RealValuedForm f,
    const std::vector<double> x_start,
    const MinimizationContext& context
);

MinimizationResult minimize_BFGS(
    RealValuedForm f,
    const std::vector<double>& x0,
    const MinimizationContext& context
);

#endif // __ANALYSIS_H__
