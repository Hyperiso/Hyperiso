#include "analysis.h"

MinimizationResult minimize_NM(RealValuedForm f, const std::vector<double> x_start, const MinimizationContext& context) {
    const std::size_t d = x_start.size();
    ScaledForm f_scaled(f, x_start);

    gsl_multimin_function gsl_f; 
    gsl_f.n = d; 
    gsl_f.f = &unwrap_lambda_multidim;
    gsl_f.params = &f_scaled;

    std::unique_ptr<gsl_multimin_fminimizer, void(*)(gsl_multimin_fminimizer*)>
        s(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, d),
                                    gsl_multimin_fminimizer_free);

    gsl_vector* x = gsl_vector_alloc(d);
    gsl_vector* step = gsl_vector_alloc(d);
    gsl_vector_set_all(x, 0.0);
    gsl_vector_set_all(step, context.step_size);

    int st = gsl_multimin_fminimizer_set(s.get(), &gsl_f, x, step);

    if (st)
        throw std::runtime_error(std::string("gsl_multimin_fminimizer_set: ") + gsl_strerror(st));

    std::size_t iter=0; int status; double size;
    do {
        ++iter; status = gsl_multimin_fminimizer_iterate(s.get());
        if (status) break;
        size = gsl_multimin_fminimizer_size(s.get());
        status = gsl_multimin_test_size(size, context.tol);
        // std::cerr << "[minimize] iter=" << iter << " status=" << gsl_strerror(status) << " size=" << size << std::endl;
    } while (status == GSL_CONTINUE && iter < context.max_iter);

    gsl_vector* gsl_argmin = gsl_multimin_fminimizer_x(s.get());
    std::vector<double> argmin(d);
    for (std::size_t i = 0; i < d; ++i) 
        argmin[i] = f_scaled.x0[i] + f_scaled.s[i] * gsl_vector_get(gsl_argmin, i);

    MinimizationResult mr;
    mr.status = status;
    mr.argmin = argmin;
    mr.min = gsl_multimin_fminimizer_minimum(s.get());

    gsl_vector_free(x);
    gsl_vector_free(step);

    return mr;
}

MinimizationResult minimize_BFGS(RealValuedForm f, const std::vector<double> &x0, const MinimizationContext &context) {
    const std::size_t d = x0.size();
    gsl_multimin_function_fdf gsl_f; 

    ScaledForm f_scaled(f, x0);

    gsl_f.n = d; 
    gsl_f.f = &unwrap_lambda_multidim;
    gsl_f.df = &unwrap_lambda_gradient;
    gsl_f.fdf = &unwrap_lambda_func_and_gradient;
    gsl_f.params = &f_scaled;

    std::unique_ptr<gsl_multimin_fdfminimizer, void(*)(gsl_multimin_fdfminimizer*)>
        s(gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, d),
                                    gsl_multimin_fdfminimizer_free);

    gsl_vector* x = gsl_vector_alloc(d);
    gsl_vector_set_all(x, 0);
    int st = gsl_multimin_fdfminimizer_set(s.get(), &gsl_f, x, context.step_size, context.line_search_tol);

    if (st)
        throw std::runtime_error(std::string("gsl_multimin_fminimizer_set: ") + gsl_strerror(st));

    std::size_t iter = 0; int status; gsl_vector* g;
    do {
        ++iter; status = gsl_multimin_fdfminimizer_iterate(s.get());
        if (status) break; // cannot improve
        g = gsl_multimin_fdfminimizer_gradient(s.get());
        status = gsl_multimin_test_gradient(g, context.tol);
    } while (status == GSL_CONTINUE && iter < context.max_iter);

    std::cout << "Minimization converged in " << iter << " iterations." << std::endl;
    gsl_vector* gsl_argmin = gsl_multimin_fdfminimizer_x(s.get());
    std::vector<double> argmin(d);
    for (std::size_t i = 0; i < d; ++i) 
        argmin[i] = f_scaled.x0[i] + f_scaled.s[i] * gsl_vector_get(gsl_argmin, i);

    MinimizationResult mr;
    mr.status = status;
    mr.argmin = argmin;
    mr.min = gsl_multimin_fdfminimizer_minimum(s.get());

    std::cout << "Minimum value = " << mr.min << std::endl;

    gsl_vector_free(x);

    return mr;
}
