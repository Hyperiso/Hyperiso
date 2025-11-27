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
    MonteCarloPredictor2(const std::shared_ptr<IModel>& model, const INuisanceSampler& sampler,
    std::map<ParamId, double> eta_mean, std::map<ParamId, std::map<ParamId, double>> eta_cov, MCPredictConfig2 cfg)
    : model_(model), sampler_(sampler), mu_(std::move(eta_mean)), Sigma_(std::move(eta_cov)), cfg_(cfg) {}
    
    MonteCarloPredictor2(const std::shared_ptr<IModel>& model, const INuisanceSampler& sampler,
    Vec eta_mean, Matrix eta_cov, MCPredictConfig2 cfg)
    : model_(model), sampler_(sampler), mu2_(std::move(eta_mean)), Sigma2_(std::move(eta_cov)), cfg_(cfg) {}

    Samples_old sample_predictions(const Vec& p, std::mt19937& rng) const {
        Samples_old out; out.reserve(cfg_.draws);
        for (std::size_t s=0; s<cfg_.draws; ++s) {
            Vec eta = sampler_.sample(mu2_, Sigma2_, rng);
            std::vector<double> value = model_->predict(p, eta);
            for (auto val : value) {
                std::cout << val << std::endl;
            }
            out.emplace_back(value);
        }
        return out;
    }

    Samples sample_predictions(const std::map<ParamId, double>& p, std::mt19937& rng) const {
        Samples out; 
        // out.reserve(cfg_.draws);
        for (std::size_t s=0; s<cfg_.draws; ++s) {
            std::cout << "s : " << s << std::endl;
            std::map<ParamId, double> eta = sampler_.sample(mu_, Sigma_, rng);
            std::map<ObservableId, double> value = model_->predict(p, eta);
            for (auto val : value) {
                std::cout << val.first.str() << " : " << val.second << std::endl;
            }
            out.emplace_back(value);
        }
        return out;
    }


    std::vector<GaussianSummary> summarize(const Vec& p, std::mt19937& rng) const {
        return gaussian_or_warn(sample_predictions(p, rng), cfg_.skew_abs_threshold);
    }

    std::vector<GaussianSummary> summarize(const std::map<ParamId, double>& p, std::mt19937& rng) const {
        std::vector<double> pa;
        for (auto& elem : p) {
            pa.push_back(elem.second);
        }
        return gaussian_or_warn(sample_predictions(p, rng), cfg_.skew_abs_threshold);
    }
private:
    const std::shared_ptr<IModel>& model_;
    const INuisanceSampler& sampler_;
    std::map<ParamId, double> mu_; std::map<ParamId, std::map<ParamId, double>> Sigma_; MCPredictConfig2 cfg_;
    Vec mu2_; Matrix Sigma2_;
};