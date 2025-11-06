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
#include "RNGHelper.h"
#include "IDecomposition.h"
#include "IDistribution.h"




class CorrelationMatrixValidator {
public:
    void validate(const Matrix& R) const;
};

// Cholesky 
class CholeskyDecomposition final : public IDecomposition {
public:
    Matrix factorize(const Matrix& R) override;
};





class RandomVectorGenerator {
public:
    RandomVectorGenerator(std::unique_ptr<IDistribution> dist,
                          std::unique_ptr<IDecomposition> decomp);

    // y = L * z, z ~ i.i.d. (E=0, Var=1). Cov(y) = L L^T = R.
    Vector generate(const Matrix& correlation) const;

private:
    std::unique_ptr<IDistribution> dist_;
    std::unique_ptr<IDecomposition> decomp_;
};



#endif