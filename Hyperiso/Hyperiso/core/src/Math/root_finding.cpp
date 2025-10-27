#include <cmath>
#include <functional>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <vector>
#include "Math.h"

bool find_bracket(const std::function<double(double)>& f,
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

double brent_root(const std::function<double(double)>& f,
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