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
//TODO : checkkkkkk
MarginalConfig MarginalConfigFactory::create(ExperimentObs oid,
                                             MarginalType marginal) {
    std::map<ExperimentObs, double> sigma;
    std::map<ExperimentObs, MarginalConfig> out;
    switch (marginal) {
    case MarginalType::GAUSSIAN:
        sigma = p(oid, DataType::STD_COMBINED);
        for (auto s : sigma) {
            out[s.first] = GaussianMarginalCfg (0.0, s.second);
        }
        return out;
        break;
    case MarginalType::FLAT:
        sigma = p(oid, DataType::STD_COMBINED);
        for (auto s : sigma) {
            out[s.first] = FlatMarginalCfg {-s.second * std::sqrt(3), s.second * std::sqrt(3)};
        }
        return out;
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
