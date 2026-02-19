#include "MarginalConfigFactory.h"

MarginalConfig MarginalConfigFactory::create(ParamId pid, MarginalType marginal) {
    double mu, sigma;

    switch (marginal) {
    case MarginalType::GAUSSIAN:
        mu = p(pid, DataType::VALUE);
        sigma = p(pid, DataType::STD_COMBINED);
        return GaussianMarginalCfg {mu, sigma};
        break;
    case MarginalType::FLAT:
        mu = p(pid, DataType::VALUE);
        sigma = p(pid, DataType::STD_COMBINED);
        return FlatMarginalCfg {mu - sigma * std::sqrt(3), mu + sigma * std::sqrt(3)};
        break;
    case MarginalType::HALF_GAUSSIAN:
        // TODO
        throw std::runtime_error("NYI");
    case MarginalType::LIKELIHOOD:
        throw std::runtime_error("NYI");
    default:
        throw std::invalid_argument("Unknown marginal type");
    }
}

MarginalConfig MarginalConfigFactory::create(BinnedObservableId oid,
                                             MarginalType marginal) {
    double sigma;

    switch (marginal) {
    case MarginalType::GAUSSIAN:
        sigma = p(oid, DataType::STD_COMBINED);
        return GaussianMarginalCfg {0.0, sigma};
        break;
    case MarginalType::FLAT:
        sigma = p(oid, DataType::STD_COMBINED);
        return FlatMarginalCfg {-sigma * std::sqrt(3), sigma * std::sqrt(3)};
        break;
    case MarginalType::HALF_GAUSSIAN:
        // TODO
        throw std::runtime_error("NYI");
    case MarginalType::LIKELIHOOD:
        throw std::runtime_error("NYI");
    default:
        throw std::invalid_argument("Unknown marginal type");
    }
}
