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

#include "Include.h"
#include "IMarginalDistribution.h"
#include "ICopula.h"
#include "Math.h"

class JointDistribution {
public:
    JointDistribution(std::vector<std::unique_ptr<IMarginalDistribution>> marginals,
                          std::unique_ptr<ICopula> copula);

    std::vector<std::vector<double>> sample(std::size_t n) const;
    std::vector<double> sample() const;
    double logpdf(std::vector<double> x) const;
    RealMatrix curvature(std::vector<double> x) const;
    std::size_t dim();
    std::vector<double> get_stds();

private:
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals_;
    std::unique_ptr<ICopula> copula_;
};



#endif