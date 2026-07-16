#ifndef DIFFCALC_H
#define DIFFCALC_H

#include <vector>
#include "functions.h"
#include "Matrix.h"

std::vector<double> gradient(const RealValuedForm& f, const std::vector<double>& x, const std::vector<double>& scales);
std::vector<double> gradient(const ScaledForm& f, const std::vector<double>& x);
RealMatrix hessian(const RealValuedForm& f, const std::vector<double>& x, const std::vector<double>& scales);
RealMatrix hessian(const ScaledForm& f, const std::vector<double>& t);
RealMatrix inverse_hessian(const RealValuedForm& f, const std::vector<double>& x, const std::vector<double>& scales);

#endif // DIFFCALC_H
