#ifndef CORRELATED_RNG_H
#define CORRELATED_RNG_H

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "RNGHelper.h"
#include "IDecomposition.h"
#include "IDistribution.h"
#include "Include.h"



class CorrelationMatrixValidator {
public:
    void validate(const Matrix& R) const;
    void validate(const std::map<ParamId, std::map<ParamId, double>>& R) const;
};

// Cholesky 
class CholeskyDecomposition final : public IDecomposition {
public:
    Matrix factorize(const Matrix& R) override;
    std::map<ParamId, std::map<ParamId, double>> factorize(const std::map<ParamId, std::map<ParamId, double>>& R);
};





class RandomVectorGenerator {
public:
    RandomVectorGenerator(std::unique_ptr<IDistribution> dist,
                          std::unique_ptr<IDecomposition> decomp);

    // y = L * z, z ~ i.i.d. (E=0, Var=1). Cov(y) = L L^T = R.
    Vector generate(const Matrix& correlation) const;
    std::map<ParamId, double> generate(const std::map<ParamId, std::map<ParamId, double>>& correlation) const;

private:
    std::unique_ptr<IDistribution> dist_;
    std::unique_ptr<IDecomposition> decomp_;
};



#endif