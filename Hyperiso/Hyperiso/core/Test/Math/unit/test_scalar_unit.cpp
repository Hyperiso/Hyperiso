#include <cassert>
#include <iostream>
#include <type_traits>
#include <complex>
#include <numbers>
#include <vector>
#include <limits>
#include "scalar.h"

static inline const complex_t& as_c(const scalar_t& z) {
    return static_cast<const complex_t&>(z);
}

static inline double cmag(const scalar_t& z) {
    return std::abs(as_c(z));
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
    std::cout << "[unit] début des tests scalar_t...\n";

    static_assert(std::is_base_of_v<complex_t, scalar_t>);
    static_assert(std::is_copy_constructible_v<scalar_t>);
    static_assert(std::is_move_constructible_v<scalar_t>);
    static_assert(std::is_copy_assignable_v<scalar_t>);
    static_assert(std::is_move_assignable_v<scalar_t>);

    {
        scalar_t z0;
        assert(dapprox(static_cast<double>(z0), 0.0));
        scalar_t z1(1.5, -2.0);
        assert(dapprox(as_c(z1).real(), 1.5) && dapprox(as_c(z1).imag(), -2.0));

        complex_t c(3.0, 4.0);
        scalar_t z2(c);      
        assert(dapprox(as_c(z2).real(), 3.0) && dapprox(as_c(z2).imag(), 4.0));

        scalar_t z3(z2);
        assert(capprox(z2, z3));

        scalar_t z4 = std::move(z3);
        assert(dapprox(as_c(z4).real(), 3.0) && dapprox(as_c(z4).imag(), 4.0));
    }

    {
        scalar_t a(2.0, 3.0), b;
        b = a;
        assert(capprox(a, b));
        scalar_t c(7.0, -1.0);
        b = std::move(c);
        assert(dapprox(as_c(b).real(), 7.0) && dapprox(as_c(b).imag(), -1.0));
    }

    {
        scalar_t z(5.0, 1e-9);
        double d = static_cast<double>(z);
        assert(dapprox(d, 5.0));
        scalar_t z2(-3.0, 1.0); 
        assert(dapprox(static_cast<double>(z2), -3.0));
    }

    {
        scalar_t z(1.0, -2.0);
        scalar_t neg = -z;
        assert(dapprox(as_c(neg).real(), -1.0) && dapprox(as_c(neg).imag(), 2.0));
    }

    {
        scalar_t a(1.0, 2.0), b(-3.0, 0.5);
        scalar_t s = a + b;
        complex_t ref = as_c(a) + as_c(b);
        assert(dapprox(as_c(s).real(), ref.real()) && dapprox(as_c(s).imag(), ref.imag()));

        scalar_t d = a - b;
        ref = as_c(a) - as_c(b);
        assert(dapprox(as_c(d).real(), ref.real()) && dapprox(as_c(d).imag(), ref.imag()));

        scalar_t m = a * b;
        ref = as_c(a) * as_c(b);
        assert(dapprox(as_c(m).real(), ref.real()) && dapprox(as_c(m).imag(), ref.imag()));

        scalar_t q = a / b;
        ref = as_c(a) / as_c(b);
        assert(dapprox(as_c(q).real(), ref.real()) && dapprox(as_c(q).imag(), ref.imag()));
    }

    {
        scalar_t z(1.0, 2.0);
        z += scalar_t(2.0, -1.0);
        assert(capprox(z, scalar_t(3.0, 1.0)));

        z -= std::complex<double>(1.0, 3.0);
        assert(capprox(z, scalar_t(2.0, -2.0)));

        z *= std::complex<double>(0.0, 1.0);
        assert(capprox(z, scalar_t(2.0 * 0.0 - (-2.0)*1.0, 2.0*1.0 + (-2.0)*0.0))); // (2,-2) * i = (2i -2i^2)?? (x+iy)*i = -y + ix -> (2,-2)*i = 2 + 2i
        assert(capprox(z, scalar_t(2.0, 2.0)));

        complex_t before_div = as_c(z);
        z /= std::complex<double>(1.0, -1.0);
        complex_t ref = before_div / complex_t(1.0, -1.0);
        assert(dapprox(as_c(z).real(), ref.real()) && dapprox(as_c(z).imag(), ref.imag()));
    }

    {
        scalar_t z(1.2, -0.5);
        assert(capprox(z + 2.0, scalar_t(3.2, -0.5)));
        assert(capprox(2.0 + z, scalar_t(3.2, -0.5)));
        assert(capprox(z - 2.0, scalar_t(-0.8, -0.5)));
        assert(capprox(2.0 - z, scalar_t(0.8, 0.5)));
        assert(capprox(z * 2.0, scalar_t(2.4, -1.0)));
        assert(capprox(2.0 * z, scalar_t(2.4, -1.0)));
        assert(capprox(z / 2.0, scalar_t(0.6, -0.25)));

        scalar_t q = 2.0 / z;
        complex_t ref = complex_t(2.0, 0.0) / as_c(z);
        assert(dapprox(as_c(q).real(), ref.real()) && dapprox(as_c(q).imag(), ref.imag()));
    }

    {
        scalar_t z(0.5, 1.0);
        std::complex<double> c(2.0, -3.0);

        scalar_t s1 = z + c;
        complex_t r1 = as_c(z) + c;
        assert(dapprox(as_c(s1).real(), r1.real()) && dapprox(as_c(s1).imag(), r1.imag()));

        scalar_t s2 = c + z;
        complex_t r2 = c + as_c(z);
        assert(dapprox(as_c(s2).real(), r2.real()) && dapprox(as_c(s2).imag(), r2.imag()));

        scalar_t p1 = z * c;
        complex_t pr1 = as_c(z) * c;
        assert(dapprox(as_c(p1).real(), pr1.real()) && dapprox(as_c(p1).imag(), pr1.imag()));

        scalar_t d1 = c / z; // instancie le template operator/(T, scalar_t)
        complex_t dr1 = c / as_c(z);
        assert(dapprox(as_c(d1).real(), dr1.real()) && dapprox(as_c(d1).imag(), dr1.imag()));
    }

    {
        scalar_t z(0.3, -0.7);
        scalar_t pcc = pow(z, scalar_t(1.2, -0.4));
        complex_t ref1 = std::pow(as_c(z), complex_t(1.2, -0.4));
        assert(dapprox(as_c(pcc).real(), ref1.real()) && dapprox(as_c(pcc).imag(), ref1.imag()));

        scalar_t pcd = pow(z, 2.5);
        complex_t ref2 = std::pow(as_c(z), 2.5);
        assert(dapprox(as_c(pcd).real(), ref2.real()) && dapprox(as_c(pcd).imag(), ref2.imag()));

        scalar_t pci = pow(z, 3);
        complex_t ref3 = std::pow(as_c(z), 3.0);
        assert(dapprox(as_c(pci).real(), ref3.real()) && dapprox(as_c(pci).imag(), ref3.imag()));
    }

    {
        scalar_t z(1.1, -0.9);

        auto check_unary = [&](auto fstd, auto fmy) {
            scalar_t r = fmy(z);
            complex_t ref = fstd(as_c(z));
            assert(dapprox(as_c(r).real(), ref.real()) && dapprox(as_c(r).imag(), ref.imag()));
        };

        check_unary([](const complex_t& c){ return std::sqrt(c); },
                    [](const scalar_t& s){ return ::sqrt(s); });
        check_unary([](const complex_t& c){ return std::sin(c); },
                    [](const scalar_t& s){ return ::sin(s); });
        check_unary([](const complex_t& c){ return std::cos(c); },
                    [](const scalar_t& s){ return ::cos(s); });
        check_unary([](const complex_t& c){ return std::tan(c); },
                    [](const scalar_t& s){ return ::tan(s); });
        check_unary([](const complex_t& c){ return std::sinh(c); },
                    [](const scalar_t& s){ return ::sinh(s); });
        check_unary([](const complex_t& c){ return std::cosh(c); },
                    [](const scalar_t& s){ return ::cosh(s); });
        check_unary([](const complex_t& c){ return std::tanh(c); },
                    [](const scalar_t& s){ return ::tanh(s); });
        check_unary([](const complex_t& c){ return std::exp(c); },
                    [](const scalar_t& s){ return ::exp(s); });
        check_unary([](const complex_t& c){ return std::log(c); },
                    [](const scalar_t& s){ return ::log(s); });
        check_unary([](const complex_t& c){ return std::asin(c); },
                    [](const scalar_t& s){ return ::asin(s); });
        check_unary([](const complex_t& c){ return std::acos(c); },
                    [](const scalar_t& s){ return ::acos(s); });
        check_unary([](const complex_t& c){ return std::atan(c); },
                    [](const scalar_t& s){ return ::atan(s); });

        {
            scalar_t a = abs(z);
            double aref = std::abs(as_c(z));
            assert(dapprox(static_cast<double>(a), aref));
            assert(dapprox(as_c(a).imag(), 0.0));

            scalar_t g = arg(z);
            double gref = std::arg(as_c(z));
            assert(dapprox(static_cast<double>(g), gref));
            assert(dapprox(as_c(g).imag(), 0.0));

            scalar_t n = norm(z);
            double nref = std::norm(as_c(z));
            assert(dapprox(static_cast<double>(n), nref));
            assert(dapprox(as_c(n).imag(), 0.0));
        }
    }

    {
        scalar_t z(-2.3, 4.5);
        assert(dapprox(real(z), -2.3));
        assert(dapprox(imag(z), 4.5));
    }

    {
        scalar_t z(0.2, 0.3);
        scalar_t r = sqrt(z);
        complex_t ref = std::sqrt(as_c(z));
        assert(dapprox(as_c(r).real(), ref.real()) && dapprox(as_c(r).imag(), ref.imag()));
    }

    std::cout << "All unit tests passed.\n";
    return 0;
}
