#include "functions.h"

ScaledForm::ScaledForm(const RealValuedForm& f, const std::vector<double>& x0, const std::vector<double>& scales) : f(f), x0(x0), s(scales) {}

double ScaledForm::operator()(const std::vector<double>& t) const {
    std::vector<double> x(t.size());
    for (size_t i = 0; i < t.size(); i++) {
        x[i] = this->x0[i] + this->s[i] * t[i];

        // std::cout << "t_i = " << t[i] << ", s_i = " << s[i] << ", x0_i = " << x0[i] << ", x_i = " << x[i] << std::endl;
    }

    double fx = this->f(x);

    // std::cout << "f(x) = " << fx << std::endl;
    return fx;
}   