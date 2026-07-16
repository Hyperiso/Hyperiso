#include "diffcalc.h"

std::vector<double> gradient(const RealValuedForm& f, const std::vector<double>& x, const std::vector<double>& scales) {
    ScaledForm f_scaled(f, x, scales);
    std::vector<double> g = gradient(f_scaled, std::vector<double>(x.size(), 0.0));
    
    for (size_t i = 0; i < g.size(); i++)
        g[i] /= f_scaled.s[i];

    return g;
}

std::vector<double> gradient(const ScaledForm &f, const std::vector<double> &t) {
    std::vector<double> grad = std::vector(t.size(), 0.0);

    double eps = std::numeric_limits<double>::epsilon();
    double h = std::cbrt(eps);
    for (size_t i = 0; i < t.size(); i++) {
        auto t_p = std::vector(t);
        auto t_m = std::vector(t);
        t_p[i] += h;
        t_m[i] -= h;
        grad[i] = (f(t_p) - f(t_m)) / (2 * h);
    }
    
    return grad;
}

RealMatrix hessian(const RealValuedForm& f, const std::vector<double>& x, const std::vector<double>& scales) {
    ScaledForm f_scaled(f, x, scales);
    RealMatrix H = hessian(f_scaled, std::vector<double>(x.size(), 0.0));
    
    for (size_t i = 0; i < H.rows(); i++) 
        for (size_t j = 0; j < H.cols(); j++)
            H.at(i, j) /= (f_scaled.s[i] * f_scaled.s[j]);

    return H;
}

RealMatrix hessian(const ScaledForm& f, const std::vector<double>& t) {
    size_t dim = t.size();
    RealMatrix H (dim, dim);

    double eps = std::numeric_limits<double>::epsilon();
    double h = std::pow(eps, 0.25);
    const double ft = f(t);
    std::vector<double> t_shift;
    std::vector<std::pair<int, int>> signs = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
    for (size_t i = 0; i < dim; i++) {
        double Hii = 0;
        for (int sign : {-1, 1}) {
            t_shift = t;
            t_shift[i] += sign * h;
            // for (double t : t_shift) std::cout << t << " ";
            // std::cout << std::endl;
            Hii += f(t_shift);
        }
        Hii = (Hii - 2 * ft) / std::pow(h, 2);

        if (!std::isfinite(Hii))
            throw std::runtime_error("Invalid value found in hessian");

        H.at(i, i) = Hii;

        for (size_t j = i + 1; j < dim; j++) {
            double Hij = 0;
            for (auto& ss : signs) {
                t_shift = t;
                t_shift[i] += ss.first * h;
                t_shift[j] += ss.second * h;
                // for (double t : t_shift) std::cout << t << " ";
                // std::cout << std::endl;
                // std::cout << f(t_shift) << std::endl;
                Hij += ss.first * ss.second * f(t_shift);
            }
            
            Hij /= 4 * h * h;

            if (!std::isfinite(Hij))
                throw std::runtime_error("Invalid value found in hessian");

            H.at(i, j) = Hij;
            H.at(j, i) = Hij;
        }
    }
    
    return 0.5 * (H + H.transpose());
}

RealMatrix inverse_hessian(const RealValuedForm& f, const std::vector<double>& x, const std::vector<double>& scales) {
    ScaledForm f_scaled(f, x, scales);
    RealMatrix H = hessian(f_scaled, std::vector<double>(x.size(), 0.0));
    RealMatrix inv_H = H.inv();

    for (size_t i = 0; i < inv_H.rows(); i++)
        for (size_t j = 0; j < inv_H.cols(); j++)
            inv_H.at(i, j) *= f_scaled.s[i] * f_scaled.s[j];
    
    return inv_H;
}