#include "DistributionFactory.h"
#include "StandardNormal.h"

std::unique_ptr<IDistribution> DistributionFactory::create(const std::string& name,
                                                unsigned int seed) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (lower == "gaussian" || lower == "normal" || lower == "gauss") {
        return std::make_unique<StandardNormal>(seed);
    }

    throw std::invalid_argument("Unkwown distribution: " + name +
                                " (try: gaussian|normal)");
}