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

    printf("Tolerance: %.2e\n", context.tol);
    printf("max_iter: %i\n", context.max_iter);

    std::unique_ptr<gsl_multimin_fminimizer, void(*)(gsl_multimin_fminimizer*)>
        s(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, d),
                                    gsl_multimin_fminimizer_free);

    if (!context.step_sizes.empty() && context.step_sizes.size() != d) {
        throw std::invalid_argument("minimize_scaled: step_sizes dimension mismatch");
    }

    gsl_vector* x = gsl_vector_alloc(d);
    gsl_vector* step = gsl_vector_alloc(d);
    for (std::size_t i = 0; i < d; ++i) { 
        gsl_vector_set(x, i, x_start[i]); 
        gsl_vector_set(step, i, context.step_sizes[i]); 

        // printf("eta_0 = %.4e, step = %.4e\n", x_start[i], context.step_sizes[i]);
    }

    int st = gsl_multimin_fminimizer_set(s.get(), &gsl_f, x, step);

    if (st) {
        throw std::runtime_error(std::string("gsl_multimin_fminimizer_set: ") + gsl_strerror(st));
    }

    std::size_t iter=0; int status; double size;
    do {
        ++iter; status = gsl_multimin_fminimizer_iterate(s.get());
        if (status) break; // cannot improve
        size = gsl_multimin_fminimizer_size(s.get());
        status = gsl_multimin_test_size(size, context.tol);
        std::cerr << "[minimize] iter=" << iter
          << " status=" << gsl_strerror(status)
          << " size=" << size
          << std::endl;
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

static inline std::vector<double> zeros(std::size_t n) {
    return std::vector<double>(n, 0.0);
}

MinimizationResult minimize_scaled(
    RealValuedForm f,
    const std::vector<double>& x0,
    const MinimizationContext& context
) {
    const std::size_t d = x0.size();

    if (!context.step_sizes.empty() && context.step_sizes.size() != d) {
        throw std::invalid_argument("minimize_scaled: step_sizes dimension mismatch");
    }

    // s_i : scale (avoid 0)
    std::vector<double> s(d, 1.0);
    if (!context.step_sizes.empty()) {
        for (std::size_t i = 0; i < d; ++i) {
            double si = context.step_sizes[i];
            if (!std::isfinite(si) || si <= 0.0) si = 1.0;
            s[i] = si;
        }
    }

    // g(u) = f(x0 + s ⊙ u)
    RealValuedForm g = [f, x0, s](std::vector<double> u) -> double {
        if (u.size() != x0.size()) {
            throw std::invalid_argument("minimize_scaled: u dimension mismatch");
        }
        std::vector<double> x(u.size());
        for (std::size_t i = 0; i < u.size(); ++i) {
            x[i] = x0[i] + s[i] * u[i];
        }
        double v = f(x);
        return std::isfinite(v) ? v : 1e300;
    };

    // Context for u u : pas = 1
    MinimizationContext uctx = context;
    uctx.step_sizes.assign(d, 1.0);

    //  u = 0 ( x = x0  u=0)
    MinimizationResult ures = minimize(g, zeros(d), uctx);

    // argmin u -> x
    std::vector<double> xhat(d);
    for (std::size_t i = 0; i < d; ++i) {
        xhat[i] = x0[i] + s[i] * ures.argmin[i];
    }

    MinimizationResult xres;
    xres.status = ures.status;
    xres.argmin = std::move(xhat);
    xres.min = f(xres.argmin); // evaluate f at real point
    
    return xres;
}