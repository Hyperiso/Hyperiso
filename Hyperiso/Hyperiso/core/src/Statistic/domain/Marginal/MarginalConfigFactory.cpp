#include "MarginalConfigFactory.h"

#include <cmath>
#include <stdexcept>
#include <utility>

MarginalConfigFactory::MarginalConfigFactory(
    std::shared_ptr<IStatParameterProxy> parameter_proxy
)
    : parameter_proxy_(std::move(parameter_proxy))
{
    if (!parameter_proxy_) {
        throw std::invalid_argument(
            "MarginalConfigFactory: parameter_proxy is null"
        );
    }
}


MarginalConfig MarginalConfigFactory::create(ParamId pid, MarginalType marginal) const {
    double mu, sigma;

    switch (marginal) {
    case MarginalType::GAUSSIAN:
        mu = (*parameter_proxy_)(pid, DataType::VALUE);
        sigma = (*parameter_proxy_)(pid, DataType::STD_COMBINED);
        return GaussianMarginalCfg {mu, sigma};
        break;
    case MarginalType::FLAT:
        mu = (*parameter_proxy_)(pid, DataType::VALUE);
        sigma = (*parameter_proxy_)(pid, DataType::STD_COMBINED);
        return FlatMarginalCfg {mu - sigma * std::sqrt(3), mu + sigma * std::sqrt(3)};
        break;
    case MarginalType::HALF_GAUSSIAN:
        // MAJ : Update data structure to store asymmetric uncertainty
        throw std::runtime_error("NYI");
    case MarginalType::LIKELIHOOD:
        // MAJ
        throw std::runtime_error("NYI");
    default:
        throw std::invalid_argument("Unknown marginal type");
    }
}


MarginalConfig MarginalConfigFactory::create(ParamId pid,
                                             MarginalType marginal,
                                             const NuisanceSpec& spec) const {
    if (spec.param_id.block != pid.block || spec.param_id.code != pid.code) {
        throw std::invalid_argument(
            "MarginalConfigFactory: nuisance specification does not match parameter block/code"
        );
    }

    if (marginal == MarginalType::FLAT) {
        const auto [lower, upper] = spec.bounds;
        if (!std::isfinite(lower) || !std::isfinite(upper) || !(lower < upper)) {
            throw std::invalid_argument(
                "MarginalConfigFactory: flat nuisance bounds must be finite and strictly ordered"
            );
        }
        return FlatMarginalCfg{lower, upper};
    }

    return create(pid, marginal);
}

MarginalConfig MarginalConfigFactory::create(ExperimentObs oid,
                                             MarginalType marginal) const {
    std::map<ExperimentObs, double> sigma =
        (*parameter_proxy_)(oid.obs, DataType::STD_COMBINED);
    MarginalConfig out;
    switch (marginal) {
    case MarginalType::GAUSSIAN:
        for (auto s : sigma) {
            if (s.first == oid) {
                out = GaussianMarginalCfg (0.0, s.second);
            }
        }
        return out;
        break;
    case MarginalType::FLAT:
        for (auto s : sigma) {
            if (s.first == oid) {
                out = FlatMarginalCfg {-s.second * std::sqrt(3), s.second * std::sqrt(3)};
            }
        }
        return out;
        break;
    case MarginalType::HALF_GAUSSIAN:
        // MAJ : Same as above
        throw std::runtime_error("NYI");
    case MarginalType::LIKELIHOOD:
        throw std::runtime_error("NYI");
    default:
        throw std::invalid_argument("Unknown marginal type");
    }
}
