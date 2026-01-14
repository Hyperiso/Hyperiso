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
#include "Include.h"
#include "IMarginalDistribution.h"
#include "ICopula.h"

class JointDistribution {
public:
    JointDistribution(std::vector<std::unique_ptr<IMarginalDistribution>> marginals,
                          std::unique_ptr<ICopula> copula);

    std::vector<Vector> sample(std::size_t n) const;
    Vector sample() const;
    double logpdf(Vector x) const;
    std::size_t dim();

private:
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals_;
    std::unique_ptr<ICopula> copula_;
};



#endif