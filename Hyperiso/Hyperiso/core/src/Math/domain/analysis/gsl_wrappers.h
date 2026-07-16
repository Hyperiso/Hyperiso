#ifndef GSL_WRAPPERS_H
#define GSL_WRAPPERS_H

#include <gsl/gsl_vector.h>
#include "functions.h"
#include "diffcalc.h"

double unwrap_lambda_unidim(double x, void *p);
double unwrap_lambda_multidim(const gsl_vector *t, void *params);
void unwrap_lambda_gradient(const gsl_vector *t, void *params, gsl_vector *g);
void unwrap_lambda_func_and_gradient(const gsl_vector *t, void *params, double *f, gsl_vector *g);

#endif // GSL_WRAPPERS_H
