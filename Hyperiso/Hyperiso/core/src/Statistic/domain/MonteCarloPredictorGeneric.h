#pragma once
#include <vector>
#include <random>
#include "INuisanceSampler.h"
#include "ports/IModel.h"
#include "Statistics.h"
#include "GaussianApprox.h"


struct MCPredictConfig2 { std::size_t draws=10000; double skew_abs_threshold=0.2; };


class MonteCarloPredictor2 {
public:
    MonteCarloPredictor2(const IModel& model, const INuisanceSampler& sampler,
    Vec eta_mean, Matrix eta_cov, MCPredictConfig2 cfg)
    : model_(model), sampler_(sampler), mu_(std::move(eta_mean)), Sigma_(std::move(eta_cov)), cfg_(cfg) {}


    Samples sample_predictions(const Vec& p, std::mt19937& rng) const {
        Samples out; out.reserve(cfg_.draws);
        for (std::size_t s=0; s<cfg_.draws; ++s) {
            Vec eta = sampler_.sample(mu_, Sigma_, rng);
            std::vector<double> value = model_.predict(p, eta);
            for (auto val : value) {
                std::cout << val << std::endl;
            }
            out.emplace_back(value);
        }
        return out;
    }


    std::vector<GaussianSummary> summarize(const Vec& p, std::mt19937& rng) const {
        return gaussian_or_warn(sample_predictions(p, rng), cfg_.skew_abs_threshold);
    }
private:
    const IModel& model_;
    const INuisanceSampler& sampler_;
    Vec mu_; Matrix Sigma_; MCPredictConfig2 cfg_;
};