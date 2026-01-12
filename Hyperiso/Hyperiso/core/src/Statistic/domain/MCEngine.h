#pragma once
#include <vector>
#include <random>
#include "INuisanceSampler.h"
#include "ports/IModel.h"
#include "Statistics.h"
#include "GaussianApprox.h"

struct MCConfig { std::size_t draws=10000; double skew_abs_threshold=0.2; };

class MonteCarloEngine {
public:
    MonteCarloEngine(const std::shared_ptr<IModel>& model, const INuisanceSampler& sampler, MCConfig cfg)
    : model_(model), sampler_(sampler), cfg_(cfg) {}

    Samples sample_predictions(const std::map<ParamId, double>& p) const {
        Samples out; 
        for (std::size_t s = 0; s < cfg_.draws; ++s) {
            std::map<ParamId, double> eta = sampler_.sample();
            std::map<ObservableId, double> value = model_->predict_optimized(p, eta);
            out.emplace_back(std::move(value));
        }
        return out;
    }

    std::vector<GaussianSummary> summarize(const std::map<ParamId, double>& p) const {
        return gaussian_fit(sample_predictions(p), cfg_.skew_abs_threshold);
    }

private:
    const std::shared_ptr<IModel>& model_;
    const INuisanceSampler& sampler_;
    MCConfig cfg_;
};


// Samples_old sample_predictions(const Vec& p, std::mt19937& rng) const {
    //     Samples_old out; out.reserve(cfg_.draws);
    //     for (std::size_t s=0; s<cfg_.draws; ++s) {
    //         Vec eta = sampler_.sample(mu2_, Sigma2_, rng);
    //         std::vector<double> value = model_->predict(p, eta);
    //         for (auto val : value) {
    //             std::cout << val << std::endl;
    //         }
    //         out.emplace_back(value);
    //     }
    //     return out;
    // }