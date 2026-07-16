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
    case MarginalType::FLAT:
        mu = (*parameter_proxy_)(pid, DataType::VALUE);
        sigma = (*parameter_proxy_)(pid, DataType::STD_COMBINED);
        return FlatMarginalCfg {mu - sigma * std::sqrt(3), mu + sigma * std::sqrt(3)};
    case MarginalType::HALF_GAUSSIAN:
        throw std::logic_error(
            "MarginalConfigFactory: HALF_GAUSSIAN parameter marginals are not implemented"
        );
    case MarginalType::LIKELIHOOD:
        throw std::logic_error(
            "MarginalConfigFactory: LIKELIHOOD parameter marginals are not implemented"
        );
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
    switch (marginal) {
    case MarginalType::GAUSSIAN:
        for (const auto& [experiment, uncertainty] : sigma) {
            if (experiment == oid) {
                return GaussianMarginalCfg {0.0, uncertainty};
            }
        }
        throw std::out_of_range(
            "MarginalConfigFactory: no combined uncertainty found for the requested observable"
        );
    case MarginalType::FLAT:
        for (const auto& [experiment, uncertainty] : sigma) {
            if (experiment == oid) {
                return FlatMarginalCfg {
                    -uncertainty * std::sqrt(3),
                    uncertainty * std::sqrt(3)
                };
            }
        }
        throw std::out_of_range(
            "MarginalConfigFactory: no combined uncertainty found for the requested observable"
        );
    case MarginalType::HALF_GAUSSIAN:
        throw std::logic_error(
            "MarginalConfigFactory: HALF_GAUSSIAN observable marginals are not implemented"
        );
    case MarginalType::LIKELIHOOD:
        throw std::logic_error(
            "MarginalConfigFactory: LIKELIHOOD observable marginals are not implemented"
        );
    default:
        throw std::invalid_argument("Unknown marginal type");
    }
}
