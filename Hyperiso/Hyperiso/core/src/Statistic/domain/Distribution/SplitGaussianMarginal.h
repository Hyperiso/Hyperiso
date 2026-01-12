#ifndef __SPLITGAUSSIANMARGINAL_H__
#define __SPLITGAUSSIANMARGINAL_H__

#include <random>
#include "IMarginalDistribution.h"
#include "Math.h"
#include "AbstractConfig.h"

struct SplitGaussianMarginalCfg : public AbstractConfig {
    double mu {0.0};
    double sigma_p {1.0};
    double sigma_m {1.0};
};

class SplitGaussianMarginal final : public IMarginalDistribution {
public:
    explicit SplitGaussianMarginal(double mu, double sigma_p, double sigma_m, unsigned int seed = std::random_device{}());

    Vector rvs(std::size_t n) override;
    double logpdf(double x) override;
    double cdf(double x) override;
    double ppf(double p) override;
    double mean() override;
    double std() override;

private:
    double mu, sigma_p, sigma_m;
    double N, w;
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};
    gsl_rng* eng_;
};

#endif // __SPLITGAUSSIANMARGINAL_H__
