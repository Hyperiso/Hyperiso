#include "Fit.h"
#include <memory>
#include <limits>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>

namespace {
struct FullPayload { const ProfiledLikelihood* like; std::size_t np; };

static double f_full(const gsl_vector* x, void* params) {
    auto* pay = static_cast<FullPayload*>(params);
    const std::size_t n  = x->size;
    const std::size_t np = pay->np;

    Vec p(np), eta(n - np);
    for (std::size_t i = 0; i < np; ++i)     p[i]   = gsl_vector_get(x, i);
    for (std::size_t j = 0; j < n - np; ++j) eta[j] = gsl_vector_get(x, np + j);

    return pay->like->ell(p, eta);
}

static double step_from(double x0) {
    // 10% relatif
    const double rel = 0.1 * std::abs(x0);

    // plancher ultra petit mais non nul (évite step=0)
    // ~ 1e-13 si x0 ~ 1e-12 ; ~ 1e-14 si x0 ~ 0 ; ~ 1e-13..1e-12 si x0 ~ 1..10
    const double floor_abs =
        1000.0 * std::numeric_limits<double>::epsilon() * (std::abs(x0) + 1.0);

    return std::max(rel, floor_abs);
}
} // namespace

FitResult MLEstimator::fit(const Vec& p0, const Vec& eta0) const {
    const std::size_t np = p0.size();
    const std::size_t ne = eta0.size();

    gsl_multimin_function f;
    f.n = np + ne;
    f.f = &f_full;

    FullPayload payload{&like_, np};
    f.params = &payload;

    std::unique_ptr<gsl_multimin_fminimizer, void(*)(gsl_multimin_fminimizer*)>
        s(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, f.n),
          gsl_multimin_fminimizer_free);

    gsl_vector* x    = gsl_vector_alloc(f.n);
    gsl_vector* step = gsl_vector_alloc(f.n);

    // init x + step
    for (std::size_t i = 0; i < np; ++i) {
        gsl_vector_set(x, i, p0[i]);
        gsl_vector_set(step, i, step_from(p0[i]));
    }
    for (std::size_t j = 0; j < ne; ++j) {
        gsl_vector_set(x, np + j, eta0[j]);
        gsl_vector_set(step, np + j, step_from(eta0[j]));
    }

    int status = gsl_multimin_fminimizer_set(s.get(), &f, x, step);
    if (status != GSL_SUCCESS) {
        std::cout << "GSL set status=" << status << " (failed to init simplex)\n";
    }

    std::size_t iter = 0;
    double size = 0.0;

    // (optionnel) logs
    // std::cout << std::setprecision(17);

    do {
        ++iter;

        status = gsl_multimin_fminimizer_iterate(s.get());
        if (status != GSL_SUCCESS) {
            std::cout << "GSL iterate status=" << status << " at iter=" << iter << "\n";
            break;
        }

        size = gsl_multimin_fminimizer_size(s.get());
        const double fval = s->fval; // valeur actuelle au "best" point
        // std::cout << "Iter " << iter << " size=" << size << " f=" << fval << "\n";

        status = gsl_multimin_test_size(size, tol);

    } while (status == GSL_CONTINUE && iter < max_iter);

    FitResult fr;
    fr.p_hat.resize(np);
    fr.eta_hat.resize(ne);

    for (std::size_t i = 0; i < np; ++i) fr.p_hat[i] = gsl_vector_get(s->x, i);
    for (std::size_t j = 0; j < ne; ++j) fr.eta_hat[j] = gsl_vector_get(s->x, np + j);

    // IMPORTANT : ell_hat cohérent avec (p_hat, eta_hat)
    fr.ell_hat = like_.ell(fr.p_hat, fr.eta_hat);

    gsl_vector_free(x);
    gsl_vector_free(step);
    return fr;
}
