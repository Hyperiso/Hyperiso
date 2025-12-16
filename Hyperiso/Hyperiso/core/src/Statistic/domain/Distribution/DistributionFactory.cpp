#include "DistributionFactory.h"
#include "StandardNormal.h"

std::unique_ptr<IDistribution> DistributionFactory::create(DistributionType name,
                                                unsigned int seed) {

    if (name == DistributionType::GAUSSIAN) {
        return std::make_unique<StandardNormal>(seed);
    }

    throw std::invalid_argument("Unkwown distribution: (try: gaussian|normal)");
}