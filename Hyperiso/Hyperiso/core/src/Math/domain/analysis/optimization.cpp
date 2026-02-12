#include "optimization.h"

bool find_bracket(const RealValuedFunction& f,
                         double x_min, double x_max,
                         double &a, double &b,
                         int n_samples)
{
    if (!(x_min < x_max)) throw std::invalid_argument("x_min < x_max required");
    std::vector<double> xs;
    xs.reserve(n_samples+1);
    for (int i = 0; i <= n_samples; ++i) {
        double t = double(i) / double(n_samples);
        xs.push_back(x_min + t * (x_max - x_min));
    }
    double fa = f(xs.front());
    if (!std::isfinite(fa)) return false;
    if (std::fabs(fa) == 0.0) { a = b = xs.front(); return true; }
    for (size_t i = 1; i < xs.size(); ++i) {
        double fb = f(xs[i]);
        if (!std::isfinite(fb)) return false;
        if (fb == 0.0) { a = b = xs[i]; return true; }
        if (fa * fb < 0.0) { a = xs[i-1]; b = xs[i]; return true; }
        fa = fb;
    }
    return false;
}

double brent_root(const RealValuedFunction& f,
                    double a, double b,
                    double xtol,
                    double ftol,
                    int max_it)
{
    double fa = f(a);
    double fb = f(b);
    if (!std::isfinite(fa) || !std::isfinite(fb))
        throw std::runtime_error("Function not finite at bracket endpoints.");
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if (fa * fb > 0.0)
        throw std::invalid_argument("Root is not bracketed (f(a) and f(b) must have opposite signs).");

    double c = a;
    double fc = fa;
    double d = b - a;
    double e = d;

    for (int iter = 0; iter < max_it; ++iter) {
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b;  fa = fb;
            b = c;  fb = fc;
            c = a;  fc = fa;
        }

        double tol = 2.0 * std::numeric_limits<double>::epsilon() * std::fabs(b) + xtol;
        double m = 0.5 * (c - b);
        if (std::fabs(m) <= tol || std::fabs(fb) <= ftol) {
            return b;
        }

        if (std::fabs(e) >= tol && std::fabs(fa) > std::fabs(fb)) {
            double s = fb / fa;
            double p, q;
            if (a == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                double r = fb / fc;
                double s2 = fa / fc;
                p = s * (2.0 * m * s2 * (s2 - r) - (b - a) * (r - 1.0));
                q = (s2 - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0) q = -q; else p = -p;
            if ( (2.0 * p) < (3.0 * m * q - std::fabs(tol * q)) && p < std::fabs(0.5 * e * q) ) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = m;
            }
        } else {
            d = m;
            e = m;
        }

        a = b; fa = fb;
        if (std::fabs(d) > tol)
            b += d;
        else
            b += (m > 0 ? tol : -tol);
        fb = f(b);
        if (!std::isfinite(fb)) throw std::runtime_error("Function returned non-finite value during iteration.");
        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a; fc = fa;
            d = b - a;
            e = d;
        }
    }

    throw std::runtime_error("brent_root: maximum iterations exceeded");
}

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