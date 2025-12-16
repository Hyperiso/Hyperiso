#include "Likelihood.h"
#include <gsl/gsl_multimin.h>
#include <memory>


namespace {
    struct ProfilePayload { const ProfiledLikelihood* self; Vec p; };

static double step_from(double x0) {
    const double rel = 0.1 * std::abs(x0);
    const double floor_abs =
        1000.0 * std::numeric_limits<double>::epsilon() * (std::abs(x0) + 1.0);
    return std::max(rel, floor_abs);
}

static double f_eta_profile(const gsl_vector* x, void* params) {
    auto* pay = static_cast<ProfilePayload*>(params);
    const std::size_t m = x->size;
    Vec eta(m);
    for (std::size_t i=0;i<m;++i) eta[i] = gsl_vector_get(x, i);
    double v = pay->self->ell(pay->p, eta);
    if (!std::isfinite(v)) return 1e300;
    return v;
}
}


double ProfiledLikelihood::ell_profiled(const Vec& p, const Vec& eta0, std::size_t max_iter, double tol) const {
    const std::size_t m = eta0.size();
    gsl_multimin_function f; f.n = m; f.f = &f_eta_profile;
    ProfilePayload payload{this, p}; f.params = &payload;


    std::unique_ptr<gsl_multimin_fminimizer, void(*)(gsl_multimin_fminimizer*)>
        s(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, m),
                                    gsl_multimin_fminimizer_free);


    gsl_vector* x = gsl_vector_alloc(m);
    gsl_vector* step = gsl_vector_alloc(m);
    for (std::size_t i=0;i<m;++i) { gsl_vector_set(x, i, eta0[i]); gsl_vector_set(step, i, step_from(eta0[i])); }


    gsl_multimin_fminimizer_set(s.get(), &f, x, step);


    std::size_t iter=0; int status; double size;
    do {
        ++iter; status = gsl_multimin_fminimizer_iterate(s.get());
        if (status) break; // cannot improve
        size = gsl_multimin_fminimizer_size(s.get());
        status = gsl_multimin_test_size(size, tol);
    } while (status == GSL_CONTINUE && iter < max_iter);


    Vec eta_hat(m);
    for (std::size_t i=0;i<m;++i) eta_hat[i] = gsl_vector_get(s->x, i);


    gsl_vector_free(x); gsl_vector_free(step);
    return ell(p, eta_hat);
}