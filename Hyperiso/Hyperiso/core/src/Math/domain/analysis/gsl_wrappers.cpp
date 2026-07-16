#include "gsl_wrappers.h"

double unwrap_lambda_unidim(double x, void *p) {
    auto fun = static_cast<std::function<double(double)>*>(p);
    return (*fun)(x);
}

double unwrap_lambda_multidim(const gsl_vector *t, void *params) {
    auto fun = static_cast<ScaledForm*>(params);
    std::vector<double> args (t->size);
    for (size_t i = 0; i < t->size; i++) {
        args[i] = gsl_vector_get(t, i);
    }
    return (*fun)(args);
}

void unwrap_lambda_gradient(const gsl_vector *t, void *params, gsl_vector *g) {
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

void unwrap_lambda_func_and_gradient(const gsl_vector *t, void *params, double *f, gsl_vector *g) {
    auto fun = static_cast<ScaledForm*>(params);
    std::vector<double> args (t->size);
    for (size_t i = 0; i < t->size; i++) {
        args[i] = gsl_vector_get(t, i);
    }
    *f = (*fun)(args);
    auto grad = gradient(*fun, args);
    for (size_t i = 0; i < grad.size(); i++) {
        // std::cout << "d_" << i << "f = " << grad[i] << std::endl;
        gsl_vector_set(g, i, grad[i]);
    }
}