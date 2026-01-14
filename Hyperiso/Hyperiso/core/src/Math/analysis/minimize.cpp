#include "analysis.h"

MinimizationResult minimize(
    RealValuedForm f,
    const std::vector<double> x_start,
    const MinimizationContext& context
) {
    const std::size_t d = x_start.size();
    gsl_multimin_function gsl_f; 
    gsl_f.n = d; 
    gsl_f.f = &unwrap_lambda_multidim;
    gsl_f.params = &f;

    std::unique_ptr<gsl_multimin_fminimizer, void(*)(gsl_multimin_fminimizer*)>
        s(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, d),
                                    gsl_multimin_fminimizer_free);

    gsl_vector* x = gsl_vector_alloc(d);
    gsl_vector* step = gsl_vector_alloc(d);
    for (std::size_t i = 0; i < d; ++i) { 
        gsl_vector_set(x, i, x_start[i]); 
        gsl_vector_set(step, i, context.step_sizes[i]); 

        printf("eta_0 = %.4e, step = %.4e\n", x_start[i], context.step_sizes[i]);
    }

    gsl_multimin_fminimizer_set(s.get(), &gsl_f, x, step);

    std::size_t iter=0; int status; double size;
    do {
        ++iter; status = gsl_multimin_fminimizer_iterate(s.get());
        if (status) break; // cannot improve
        size = gsl_multimin_fminimizer_size(s.get());
        status = gsl_multimin_test_size(size, context.tol);
    } while (status == GSL_CONTINUE && iter < context.max_iter);

    std::vector<double> argmin(d);
    for (std::size_t i = 0; i < d; ++i) 
        argmin[i] = gsl_vector_get(s->x, i);

    gsl_vector_free(x);
    gsl_vector_free(step);

    MinimizationResult mr;
    mr.status = status;
    mr.argmin = argmin;
    mr.min = f(argmin);

    return mr;
}