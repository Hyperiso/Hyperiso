#include "analysis.h"

std::vector<double> gradient(RealValuedForm f, const std::vector<double>& x) {
    std::vector<double> grad = std::vector(x.size(), 0.0);

    double h = 1e-8;
    for (size_t i = 0; i < x.size(); i++) {
        auto x_p = std::vector(x);
        auto x_m = std::vector(x);
        x_p[i] *= (1 + h);
        x_m[i] *= (1 - h);
        grad[i] = (f(x_p) - f(x_m)) / (2 * h * x[i]);
    }
    
    return grad;
}

RealMatrix hessian(RealValuedForm f, const std::vector<double>& x) {
    size_t dim = x.size();
    RealMatrix H (dim, dim);

    double h = 1e-5;
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            auto x_pp = std::vector(x);
            auto x_pm = std::vector(x);
            auto x_mp = std::vector(x);
            auto x_mm = std::vector(x);
            x_pp[i] *= (1 + h);
            x_pm[i] *= (1 + h);
            x_mp[i] *= (1 - h);
            x_mm[i] *= (1 - h);
            x_pp[j] *= (1 + h);
            x_pm[j] *= (1 - h);
            x_mp[j] *= (1 + h);
            x_mm[j] *= (1 - h);
            H.at(i, j) = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4 * h * h * x[i] * x[j]);
        }
    }
    
    return H;
}