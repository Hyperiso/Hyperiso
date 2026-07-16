// test_scalar_integration.cpp
#include <cassert>
#include <iostream>
#include <numbers>
#include <vector>
#include <complex>
#include "scalar.h"

static inline const complex_t& as_c(const scalar_t& z) {
    return static_cast<const complex_t&>(z);
}
static inline bool capprox(const scalar_t& a, const scalar_t& b, double eps=1e-12) {
    const double diff = std::abs(as_c(a) - as_c(b));
    const double scale = 1.0 + std::abs(as_c(a)) + std::abs(as_c(b));
    return diff <= eps * scale;
}
static inline bool dapprox(double a, double b, double eps=1e-12) {
    return std::abs(a - b) <= eps * (1.0 + std::abs(a) + std::abs(b));
}

int main() {
    std::cout << "[integration] début des tests...\n";

    // 1) Euler : exp(i*pi) + 1 -> 0
    {
        const double pi = std::numbers::pi;
        scalar_t i(0.0, 1.0);
        scalar_t r = exp(i * pi) + 1.0;
        assert(capprox(r, scalar_t(0.0, 0.0), 1e-12));
    }

    // 2) Quadratique : z^2 + 2z + 5 = 0  => racines -1 ± 2i
    {
        scalar_t a(1.0, 0.0), b(2.0, 0.0), c(5.0, 0.0);
        scalar_t disc = b*b - 4.0*a*c;                // Δ = 4 - 20 = -16
        scalar_t sqrt_disc = sqrt(disc);               // 4i
        scalar_t two_a = 2.0 * a;

        scalar_t z1 = (-b + sqrt_disc) / two_a;
        scalar_t z2 = (-b - sqrt_disc) / two_a;

        assert(capprox(z1, scalar_t(-1.0, 2.0)));
        assert(capprox(z2, scalar_t(-1.0, -2.0)));

        // Vérif inverse : remonter dans le polynôme
        auto P = [&](const scalar_t& z) { return z*z + b*z + c; };
        assert(capprox(P(z1), scalar_t(0.0, 0.0)));
        assert(capprox(P(z2), scalar_t(0.0, 0.0)));
    }

    // 3) Identité trig en complexe + mélange de types
    {
        scalar_t z(0.7, -1.1);
        scalar_t lhs = pow(sin(z), 2) + pow(cos(z), 2); // sin^2 + cos^2
        assert(capprox(lhs, scalar_t(1.0, 0.0)));

        // Expression mixte : (2 + i) * z  +  3  -  (z / (1 - i))
        std::complex<double> ci(2.0, 1.0);
        scalar_t expr = ci * z + 3.0 - z / std::complex<double>(1.0, -1.0);

        // Vérif num. avec complex "pur"
        complex_t ref = (complex_t(2.0,1.0) * as_c(z)) + complex_t(3.0,0.0) - (as_c(z) / complex_t(1.0,-1.0));
        assert(dapprox(as_c(expr).real(), ref.real()) && dapprox(as_c(expr).imag(), ref.imag()));
    }

    std::cout << "All integration tests passed.\n";
    return 0;
}
