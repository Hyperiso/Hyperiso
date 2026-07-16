#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <gsl/gsl_integration.h>
#include "functions.h"
#include "gsl_wrappers.h"

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

#endif // INTEGRATION_H
