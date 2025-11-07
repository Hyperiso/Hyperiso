#pragma once
#include <vector>
#include <random>
#include <stdexcept>
#include <cmath>


#include "INuisanceSampler.h"
#include "RNGHelper.h"
#include "RandomVectorGenerator.h"


// Draws y ~ Corr(R) with unit variances using your RVG, then scales: η = μ + D y,
// where Σ = D R D and D = diag(σ).
class RvgNuisanceSampler final : public INuisanceSampler {
public:
    explicit RvgNuisanceSampler(RandomVectorGenerator& rvg) : rvg_(rvg) {}

    Vec sample(const Vec& mean, const Matrix& cov, std::mt19937&) const override {
        const std::size_t n = mean.size();
        if (cov.size()!=n) throw std::invalid_argument("cov dimension mismatch");
        for (const auto& r: cov) if (r.size()!=n) throw std::invalid_argument("cov not square");


        // Build stddevs and correlation matrix
        std::vector<double> sigma(n);
        Matrix R(n, Vec(n, 0.0));
        for (std::size_t i=0;i<n;++i) {
            const double v = cov[i][i];
            if (v<=0.0) throw std::invalid_argument("non-positive variance on diag(cov)");
                sigma[i] = std::sqrt(v);
        }
        for (std::size_t i=0;i<n;++i) {
            for (std::size_t j=0;j<n;++j) {
                R[i][j] = cov[i][j] / (sigma[i]*sigma[j]);
            }
        }
        // y ~ Corr(R) with unit variances
        Vector y = rvg_.generate(R);
        // scale to covariance Σ and add mean
        Vec eta(n);
        for (std::size_t i=0;i<n;++i) eta[i] = mean[i] + sigma[i]*y[i];
        return eta;
    }
private:
    RandomVectorGenerator& rvg_;
};