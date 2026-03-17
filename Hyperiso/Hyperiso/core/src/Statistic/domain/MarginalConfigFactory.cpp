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
std::vector<MarginalConfig> MarginalConfigFactory::create(BinnedObservableId oid,
                                             MarginalType marginal) {
    std::vector<double> sigma;
    std::vector<MarginalConfig> out;
    switch (marginal) {
    case MarginalType::GAUSSIAN:
        sigma = p(oid, DataType::STD_COMBINED);
        for (auto s : sigma) {
            out.push_back(GaussianMarginalCfg (0.0, s));
        }
        return out;
        break;
    case MarginalType::FLAT:
        sigma = p(oid, DataType::STD_COMBINED);
        for (auto s : sigma) {
            out.push_back(FlatMarginalCfg {-s * std::sqrt(3), s * std::sqrt(3)});
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
