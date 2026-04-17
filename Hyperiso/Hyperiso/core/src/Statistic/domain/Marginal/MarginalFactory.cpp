#include "MarginalFactory.h"

std::unique_ptr<IMarginalDistribution> make(unsigned int seed, MarginalConfig cfg) {
    return std::visit([seed](auto&& c) -> std::unique_ptr<IMarginalDistribution> {
        using T = std::decay_t<decltype(c)>;
        if constexpr (std::is_same_v<T, FlatMarginalCfg>)
            return std::make_unique<FlatMarginal>(c.a, c.b, seed);
        else if constexpr (std::is_same_v<T, GaussianMarginalCfg>)
            return std::make_unique<GaussianMarginal>(c.mu, c.sigma, seed);
        else if constexpr (std::is_same_v<T, SplitGaussianMarginalCfg>)
            return std::make_unique<SplitGaussianMarginal>(c.mu, c.sigma_p, c.sigma_m, seed);
        else if constexpr (std::is_same_v<T, LikelihoodMarginalCfg>)
            return std::make_unique<LikelihoodMarginal>(c.values, c.weights, seed, true);
    }, cfg);
}

std::unique_ptr<IMarginalDistribution> MarginalFactory::create(MarginalType name, MarginalConfig cfg, unsigned int seed) {
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