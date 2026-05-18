#include "ContourEngine.h"

namespace {

std::size_t count_total_points(const Contour& c) {
    std::size_t n = 0;
    for (const auto& path : c.paths) {
        n += path.size();
    }
    return n;
}

}

ContourEngine::ContourEngine(std::shared_ptr<ILikelihood> base, const ContourConfig &cfg) : cfg(cfg) {
    // TODO : Maybe allow to change fit backend, for now only Minuit hardcoded.
    // std::shared_ptr<Profiler> profiler = std::make_shared<Profiler>(fit_app::make_minuit_backend());
    
    ProfilerMode profiler_mode = cfg.profile_backend;
    std::shared_ptr<Profiler> profiler =
        std::make_shared<Profiler>(
            fit_app::make_minuit_backend(),
            profiler_mode
        );

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

// Contour ContourEngine::compute_contour(double z, std::array<double, 4> bounds, std::size_t resolution) {
//     ScalarField2D field = [this] (double x, double y) {
//         return this->likelihood.profiled_nll(x, y) - this->cfg.fr.ell_hat;
//     };

//     ContourRequest cr;
//     cr.bounds = bounds;
//     cr.level = z * z / 2; 
//     cr.resolution = resolution;

//     return this->extractor->extract(field, cr);
// }

Contour ContourEngine::compute_contour(double z, std::array<double, 4> bounds, std::size_t resolution) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();

    // ScalarField2D field = [this] (double x, double y) {
    //     return this->likelihood.profiled_nll(x, y) - this->cfg.fr.ell_hat;
    // };

    const double x_ref = cfg.fr.p_hat.at(cfg.x_id);
    const double y_ref = cfg.fr.p_hat.at(cfg.y_id);


    const double reference_nll =
        this->likelihood.profiled_nll(x_ref, y_ref);

    // ScalarField2D field = [this, reference_nll](double x, double y) {
    //     const double v = this->likelihood.profiled_nll(x, y);

    //     if (!std::isfinite(v)) {
    //         return 1e50;
    //     }

    //     return std::clamp(v - reference_nll, -1e50, 1e50);
    // };

    auto cache = std::make_shared<std::map<Point, double>>();

    ScalarField2D field = [this, reference_nll, cache](double x, double y) {
        Point key{x, y};

        auto it = cache->find(key);
        if (it != cache->end()) {
            return it->second;
        }

        const double raw = this->likelihood.profiled_nll(x, y);
        const double val = std::isfinite(raw)
            ? std::clamp(raw - reference_nll, -1e50, 1e50)
            : 1e50;

        (*cache)[key] = val;
        return val;
    };

    ContourRequest cr;
    cr.bounds = bounds;
    cr.level = z * z / 2.0;
    cr.resolution = resolution;

    auto defs = this->likelihood.get_param_defs();

    defs[0].value = cfg.fr.p_hat.at(cfg.x_id);
    defs[0].step_hint = std::max(cfg.fr.p_hat_std.at(cfg.x_id), 1e-3);
    defs[0].limits = std::make_pair(bounds[0], bounds[1]);

    defs[1].value = cfg.fr.p_hat.at(cfg.y_id);
    defs[1].step_hint = std::max(cfg.fr.p_hat_std.at(cfg.y_id), 1e-3);
    defs[1].limits = std::make_pair(bounds[2], bounds[3]);

    cr.p_defs[0] = defs[0];
    cr.p_defs[1] = defs[1];

    if (cfg.on_progress) {
        ContourProgressEvent ev;
        ev.type = ContourProgressEventType::Started;
        ev.level = cr.level;
        ev.message = "Contour computation started";
        cfg.on_progress(ev);
    }

    try {
        Contour c = this->extractor->extract(field, cr);

        if (cfg.on_progress) {
            ContourProgressEvent ev;
            ev.type = ContourProgressEventType::Finished;
            ev.level = c.level;
            ev.n_paths = c.paths.size();
            ev.n_points = count_total_points(c);
            ev.elapsed_seconds =
                std::chrono::duration<double>(clock::now() - t0).count();
            ev.message = c.success ? "Contour computation finished"
                                   : "Contour computation finished but contour is marked unsuccessful";
            cfg.on_progress(ev);
        }

        return c;
    }
    catch (const std::exception& e) {
        if (cfg.on_progress) {
            ContourProgressEvent ev;
            ev.type = ContourProgressEventType::Failed;
            ev.level = cr.level;
            ev.elapsed_seconds =
                std::chrono::duration<double>(clock::now() - t0).count();
            ev.message = e.what();
            cfg.on_progress(ev);
        }
        throw;
    }
}

// Contour ContourEngine::compute_contour(double z, std::array<double, 4> bounds, std::size_t resolution) {
//     ScalarField2D field = [this] (double x, double y) {
//         return this->likelihood.profiled_nll(x, y) - this->cfg.fr.ell_hat;
//     };

//     ContourRequest cr;
//     cr.bounds = bounds;
//     cr.level = z * z / 2.0;
//     cr.resolution = resolution;

//     //TODO : Niels -> L'erreur c'était qu'on ne def pas les p ici
//     auto defs = this->likelihood.get_param_defs();

//     defs[0].value = cfg.fr.p_hat.at(cfg.x_id);
//     defs[0].step_hint = std::max(cfg.fr.p_hat_std.at(cfg.x_id), 1e-3);
//     defs[0].limits = std::make_pair(bounds[0], bounds[1]);

//     defs[1].value = cfg.fr.p_hat.at(cfg.y_id);
//     defs[1].step_hint = std::max(cfg.fr.p_hat_std.at(cfg.y_id), 1e-3);
//     defs[1].limits = std::make_pair(bounds[2], bounds[3]);

//     cr.p_defs[0] = defs[0];
//     cr.p_defs[1] = defs[1];

//     return this->extractor->extract(field, cr);
// }

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

// std::shared_ptr<IContourExtractor> ContourEngine::build_contour_extractor(
//     ContourAlgorithm ca)
// {
//     switch (ca) {
//     case ContourAlgorithm::AMS:
//         return std::make_shared<AMSContourExtractor>();
//     case ContourAlgorithm::MINUIT:
//         return std::make_shared<MnContourExtractor>();
//     default:
//         throw std::invalid_argument("Unknown contouring algorithm.");
//     }
// }

std::shared_ptr<IContourExtractor> ContourEngine::build_contour_extractor(
    ContourAlgorithm ca)
{
    switch (ca) {
    case ContourAlgorithm::AMS:
        return std::make_shared<AMSContourExtractor>(cfg.on_progress);
    case ContourAlgorithm::MINUIT:
        return std::make_shared<MnContourExtractor>();
    default:
        throw std::invalid_argument("Unknown contouring algorithm.");
    }
}