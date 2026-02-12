#include "functions.h"

ScaledForm::ScaledForm(const RealValuedForm& f, const std::vector<double>& x0) : f(f), x0(x0) {
    this->s = std::vector(x0);
    std::transform(this->s.begin(), this->s.end(), this->s.begin(), std::abs<double>);
    double max_abs_x = *std::max_element(this->s.begin(), this->s.end());
    double eps = std::numeric_limits<double>::epsilon();
    double sigma_min = std::sqrt(eps) * (1 + max_abs_x);
    std::transform(this->s.begin(), this->s.end(), this->s.begin(), [sigma_min] (double x) { return std::max(x, sigma_min); });
}

double ScaledForm::operator()(const std::vector<double>& t) const {
    std::vector<double> x(t.size());
    for (size_t i = 0; i < t.size(); i++) {
        x[i] = this->x0[i] + this->s[i] * t[i];

        std::cout << "t_i = " << t[i] << ", s_i = " << s[i] << ", x0_i = " << x0[i] << ", x_i = " << x[i] << std::endl;
    }

    return this->f(x);
}   