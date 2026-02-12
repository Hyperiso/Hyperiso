#ifndef __DIFFCALC_H__
#define __DIFFCALC_H__

#include <vector>
#include "functions.h"
#include "../linalg/Matrix.h"

std::vector<double> gradient(const RealValuedForm& f, const std::vector<double>& x);
std::vector<double> gradient(const ScaledForm& f, const std::vector<double>& x);
RealMatrix hessian(const RealValuedForm& f, const std::vector<double>& x);
RealMatrix hessian(const ScaledForm& f, const std::vector<double>& t);
RealMatrix inverse_hessian(const RealValuedForm& f, const std::vector<double>& x);

#endif // __DIFFCALC_H__
