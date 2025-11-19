#pragma once
#include <vector>
#include <random>
#include "LinearAlgebra.h"
#include "Statistics.h"
#include "GaussianApprox.h"
#include "ports/IModel.h"


struct MCPredictConfig {
    std::size_t draws = 10000;
    double skew_abs_threshold = 0.2; // decide Gaussian approx validity
};


class MonteCarloPredictor {
public:
    MonteCarloPredictor(IModel& model, Vec eta_mean, Matrix eta_cov, MCPredictConfig cfg)
    : model_(model), mu_(std::move(eta_mean)), Sigma_(std::move(eta_cov)), cfg_(cfg), Sigma_chol_(SPDMatrix::cholesky(Sigma_)) {}


    // Sample nuisances and return samples of O_th
    Samples sample_predictions(const Vec& p, std::mt19937& rng) {
        std::normal_distribution<double> N01(0.0, 1.0);
        const std::size_t ne = mu_.size();
        const std::size_t N = cfg_.draws;
        Samples out; out.reserve(N);
        for (std::size_t s=0;s<N;++s) {
        Vec z(ne); for (std::size_t i=0;i<ne;++i) z[i] = N01(rng);
        // eta = mu + L z
        Vec eta(ne, 0.0);
        for (std::size_t i=0;i<ne;++i) {
        double acc = 0.0; for (std::size_t k=0;k<=i;++k) acc += Sigma_chol_.L[i][k] * z[k];
        eta[i] = mu_[i] + acc;
        }
        out.emplace_back(model_.predict(p, eta));
        }
        return out;
    }


    // Return Gaussian summaries per observable and a flag for approx validity
    std::vector<GaussianSummary> summarize(const Vec& p, std::mt19937& rng) {
    auto S = sample_predictions(p, rng);
    return gaussian_or_warn(S, cfg_.skew_abs_threshold);
    }


private:
    IModel& model_;
    Vec mu_;
    Matrix Sigma_;
    MCPredictConfig cfg_;
    SPDMatrix Sigma_chol_;
};