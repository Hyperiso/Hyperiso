#include "MarginalFactory.h"

std::unique_ptr<IMarginalDistribution> make(unsigned int seed, Config cfg) {
    return std::visit([](auto&& c) -> std::unique_ptr<IMarginalDistribution> {
        using T = std::decay_t<decltype(c)>;
        if constexpr (std::is_same_v<T, FlatMarginalCfg>)
            return std::make_unique<FlatMarginal>(seed, c.a, c.b);
        else if constexpr (std::is_same_v<T, GaussianMarginalCfg>)
            return std::make_unique<GaussianMarginal>(seed, c.mu, c.sigma);
        else if constexpr (std::is_same_v<T, SplitGaussianMarginalCfg>)
            return std::make_unique<SplitGaussianMarginal>(seed, c.mu, c.sigma_p, c.sigma_m);
        else if constexpr (std::is_same_v<T, LikelihoodMarginalCfg>)
            return std::make_unique<LikelihoodMarginal>(seed, c.values, c.weights, true);
    }, cfg);
}

std::unique_ptr<IMarginalDistribution> DistributionFactory::create(MarginalType name, Config cfg, unsigned int seed) {
    switch (name) {
        case MarginalType::GAUSSIAN:
            return make(seed, cfg);
        case MarginalType::HALF_GAUSSIAN:
            return make(seed, cfg);
        case MarginalType::FLAT:
            return make(seed, cfg);
        //TODO : 
        case MarginalType::LIKELIHOOD:
            throw std::invalid_argument("LIKELIHOOD needs data (values, weights). Use createLikelihood(...).");
            return make(seed, cfg);
        default:
            throw std::invalid_argument("Unknown distribution.");
    }
}