#include "DistributionFactory.h"
#include "StandardNormal.h"
#include "StandardFlat.h"
#include "LikelihoodDiscrete.h"

std::unique_ptr<IDistribution> DistributionFactory::create(DistributionType name,
                                                unsigned int seed) {

    switch (name) {
        case DistributionType::GAUSSIAN:
            return std::make_unique<StandardNormal>(seed);
        case DistributionType::FLAT:
            return std::make_unique<StandardFlat>(seed);
        //TODO : 
        case DistributionType::LIKELIHOOD:
            throw std::invalid_argument("LIKELIHOOD needs data (values, weights). Use createLikelihood(...).");
        default:
            throw std::invalid_argument("Unknown distribution.");
    }
}