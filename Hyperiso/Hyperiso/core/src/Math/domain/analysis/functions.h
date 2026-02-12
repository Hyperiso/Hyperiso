#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <functional>
#include <limits>
#include <algorithm>
#include "scalar.h"

using RealValuedFunction = std::function<double(double)>;
using RealValuedForm = std::function<double(std::vector<double>)>;
using ComplexValuedFunction = std::function<scalar_t(double)>;

struct ScaledForm {
    RealValuedForm f;
    std::vector<double> x0;
    std::vector<double> s;

    ScaledForm(const RealValuedForm& f, const std::vector<double>& x0);
    double operator()(const std::vector<double>& t) const;
};

#endif // FUNCTIONS_H
