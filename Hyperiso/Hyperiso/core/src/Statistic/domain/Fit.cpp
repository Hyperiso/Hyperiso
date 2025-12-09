#include "Fit.h"
#include <memory>
#include <limits>


namespace {
struct FullPayload { const ProfiledLikelihood* like; std::size_t np; };


static double f_full(const gsl_vector* x, void* params) {
    auto* pay = static_cast<FullPayload*>(params);
    const std::size_t n = x->size;
    const std::size_t np = pay->np;
    Vec p(np), eta(n-np);
    for (std::size_t i=0;i<np;++i) p[i] = gsl_vector_get(x, i);
    for (std::size_t j=0;j<n-np;++j) eta[j] = gsl_vector_get(x, np+j);
    return pay->like->ell(p, eta);
}
}


FitResult MLEstimator::fit(const Vec& p0, const Vec& eta0) const {
    const std::size_t np = p0.size();
    const std::size_t ne = eta0.size();
    gsl_multimin_function f; f.n = np+ne; f.f = &f_full;
    FullPayload payload{&like_, np}; f.params = &payload;


    std::unique_ptr<gsl_multimin_fminimizer, void(*)(gsl_multimin_fminimizer*)>
        s(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, f.n),
            gsl_multimin_fminimizer_free);


    gsl_vector* x = gsl_vector_alloc(f.n);
    gsl_vector* step = gsl_vector_alloc(f.n);
    for (std::size_t i=0;i<np;++i) { gsl_vector_set(x, i, p0[i]); gsl_vector_set(step, i, 0.2); }
    for (std::size_t j=0;j<ne;++j) { gsl_vector_set(x, np+j, eta0[j]); gsl_vector_set(step, np+j, 0.2); }


    gsl_multimin_fminimizer_set(s.get(), &f, x, step);


    std::size_t iter=0; int status; double size;
    do {
        ++iter; status = gsl_multimin_fminimizer_iterate(s.get());
        if (status) break;
        size = gsl_multimin_fminimizer_size(s.get());
        status = gsl_multimin_test_size(size, tol);
    } while (status == GSL_CONTINUE && iter < max_iter);


    FitResult fr; fr.p_hat.resize(np); fr.eta_hat.resize(ne);
    for (std::size_t i=0;i<np;++i) fr.p_hat[i] = gsl_vector_get(s->x, i);
    for (std::size_t j=0;j<ne;++j) fr.eta_hat[j] = gsl_vector_get(s->x, np+j);
    fr.ell_hat = gsl_multimin_fminimizer_minimum(s.get());


    gsl_vector_free(x); gsl_vector_free(step);
    return fr;
}


std::pair<double,double> MLEstimator::confidence_interval_1d(
    const FitResult& fr,
    double p_min, double p_max,
    int grid_points,
    double alpha,
    const Vec& eta_init) const {
    const double threshold = gsl_cdf_chisq_Pinv(1.0 - alpha, 1); // dof=1 for profile of 1 param


    double left = std::numeric_limits<double>::quiet_NaN();
    double right = std::numeric_limits<double>::quiet_NaN();


    double prev_p = p_min;
    double prev_T = test_statistic(Vec{prev_p}, fr, eta_init);


    for (int i=1;i<=grid_points;++i) {
        const double p = p_min + (p_max - p_min)*i/static_cast<double>(grid_points);
        const double T = test_statistic(Vec{p}, fr, eta_init);
        // detect crossings of threshold
        if (std::isnan(left) && ((prev_T - threshold) * (T - threshold) <= 0.0)) left = prev_p;
        if ( ((prev_T - threshold) * (T - threshold) <= 0.0)) right = p;
        prev_p = p; prev_T = T;
    }


    return {left, right};
}