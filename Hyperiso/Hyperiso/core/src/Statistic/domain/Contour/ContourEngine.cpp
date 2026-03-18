#include "ContourEngine.h"

ContourEngine::ContourEngine(std::shared_ptr<ILikelihood> base, const ContourConfig &cfg) : cfg(cfg) {
    // TODO : Maybe allow to change fit backend, for now only Minuit hardcoded.
    std::shared_ptr<Profiler> profiler = std::make_shared<Profiler>(fit_app::make_minuit_backend());
    
    std::shared_ptr<IProfilingStrategy> strategy;
    if (cfg.profiling_method == ProfilingMethod::SLICE) {
        strategy = std::make_shared<SliceProfilingStrategy>(cfg.x_id, cfg.y_id, cfg.fr);
    } else {
        strategy = std::make_shared<ProjectionProfilingStrategy>(cfg.x_id, cfg.y_id, cfg.fr);
    }

    std::shared_ptr<ILikelihood> maybe_constrained_base = base;
    if (cfg.profiling_method == ProfilingMethod::PRIOR_CONSTRAINED_PROJECTION) {
        std::vector<std::size_t> constrained_idx;
        for (size_t i = 0; i < cfg.fr.p_hat.size(); i++) {
            if (i != cfg.x_id && i != cfg.y_id)
                constrained_idx.emplace_back(i);
        }
        
        maybe_constrained_base = std::make_shared<WithGaussianConstraints>(
            base, 
            this->build_constraints_distribution(),
            constrained_idx
        );
    }

    this->likelihood = ProfiledLikelihood2D(maybe_constrained_base, profiler, strategy);

    if (cfg.fallback_contour_method.has_value()) {
        this->extractor = std::make_shared<WithFallback>(
            this->build_contour_extractor(cfg.primary_contour_method),
            this->build_contour_extractor(cfg.fallback_contour_method.value())
        );
    } else {
        this->extractor = this->build_contour_extractor(cfg.primary_contour_method);
    }
}

Contour ContourEngine::compute_contour(double z, std::array<double, 4> bounds, std::size_t resolution) {
    ScalarField2D field = [this] (double x, double y) {
        return this->likelihood.profiled_nll(x, y) - this->cfg.fr.ell_hat;
    };

    ContourRequest cr;
    cr.bounds = bounds;
    cr.level = z * z / 2; // TODO : check !
    cr.resolution = resolution;

    return this->extractor->extract(field, cr);
}

std::shared_ptr<JointDistribution> ContourEngine::build_constraints_distribution() {
    std::vector<std::unique_ptr<IMarginalDistribution>> fitted_marginals;

    for (size_t i = 0; i < cfg.fr.p_hat.size(); i++) {
        if (i != cfg.x_id && i != cfg.y_id)
            fitted_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>(cfg.fr.p_hat[i], cfg.fr.p_hat_std[i])));
    }

    RealMatrix R_fitted = cfg.fr.p_hat_correlations;
    R_fitted.remove_row_and_column(cfg.y_id);
    R_fitted.remove_row_and_column(cfg.x_id);
    std::unique_ptr<ICopula> fitted_copula = std::make_unique<GaussianCopula>(std::random_device{}(), R_fitted);

    return std::make_unique<JointDistribution>(std::move(fitted_marginals), std::move(fitted_copula));
}

std::shared_ptr<IContourExtractor> ContourEngine::build_contour_extractor(
    ContourAlgorithm ca)
{
    switch (ca) {
    case ContourAlgorithm::AMS:
        return std::make_shared<AMSContourExtractor>();
    case ContourAlgorithm::MINUIT:
        return std::make_shared<MnContourExtractor>();
    default:
        throw std::invalid_argument("Unknown contouring algorithm.");
    }
}
