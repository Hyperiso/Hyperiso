#pragma once
#include <vector>
#include <random>
#include "INuisanceSampler.h"
#include "ports/IModel.h"
#include "Statistics.h"
#include "GaussianApprox.h"
#include "Indexing.h"

struct MCConfig {
    std::size_t draws = 10000;
    double skew_abs_threshold = 0.2;

    double covariance_ridge_rel = 1e-8;
    double covariance_ridge_abs = 1e-12;
};

struct MCRealization {
    ObsSamples sampled_obss;
    NuisanceSamples sampled_params;
};

struct MCObservableCovariance {
    std::vector<BinnedObservableId> ids;
    std::vector<double> mean;
    RealMatrix covariance;
    RealMatrix covariance_inv;
};



inline MCObservableCovariance covariance_from_obs_samples(
    const ObsSamples& S,
    const std::vector<BinnedObservableId>& ids,
    double ridge_rel = 1e-8,
    double ridge_abs = 1e-12
) {
    if (S.empty()) {
        throw std::invalid_argument("covariance_from_obs_samples: no samples");
    }
    if (ids.empty()) {
        throw std::invalid_argument("covariance_from_obs_samples: no observable ids");
    }
    if (S.size() < 2) {
        throw std::invalid_argument("covariance_from_obs_samples: need at least two samples");
    }

    const std::size_t N = S.size();
    const std::size_t D = ids.size();

    std::vector<double> mean(D, 0.0);
    for (const auto& row : S) {
        for (std::size_t d = 0; d < D; ++d) {
            mean[d] += row.at(ids[d]);
        }
    }
    for (double& v : mean) {
        v /= static_cast<double>(N);
    }

    RealMatrix cov(D, D);
    for (std::size_t i = 0; i < D; ++i) {
        for (std::size_t j = 0; j < D; ++j) {
            double s = 0.0;
            for (const auto& row : S) {
                s += (row.at(ids[i]) - mean[i]) * (row.at(ids[j]) - mean[j]);
            }
            cov.at(i, j) = s / static_cast<double>(N - 1);
        }
    }

    // Force symmetry.
    for (std::size_t i = 0; i < D; ++i) {
        for (std::size_t j = i + 1; j < D; ++j) {
            const double v = 0.5 * (cov.at(i, j) + cov.at(j, i));
            cov.at(i, j) = v;
            cov.at(j, i) = v;
        }
    }

    double trace = 0.0;
    for (std::size_t i = 0; i < D; ++i) {
        trace += std::max(0.0, cov.at(i, i));
    }

    const double scale = D > 0 ? trace / static_cast<double>(D) : 1.0;
    const double ridge = std::max(ridge_abs, ridge_rel * std::max(scale, 1.0));

    for (std::size_t i = 0; i < D; ++i) {
        cov.at(i, i) += ridge;
    }

    MCObservableCovariance out;
    out.ids = ids;
    out.mean = mean;
    out.covariance = cov;
    out.covariance_inv = cov.inv();
    return out;
}

inline std::vector<BinnedObservableId> covariance_ids_from_first_sample(
    const ObsSamples& S
) {
    if (S.empty()) {
        throw std::invalid_argument("covariance_ids_from_first_sample: no samples");
    }

    std::vector<BinnedObservableId> ids;
    ids.reserve(S.front().size());
    for (const auto& [id, _] : S.front()) {
        ids.push_back(id);
    }
    return ids;
}

struct MCResult {
    MCRealization mc_real;
    std::vector<GaussianSummary> summary;

    MCObservableCovariance covariance;
};

class MonteCarloEngine {
public:
    MonteCarloEngine(const std::shared_ptr<IModel>& model, const INuisanceSampler& sampler, MCConfig cfg)
    : model_(model), sampler_(sampler), cfg_(cfg) {}

    MCRealization sample_predictions(const std::map<ParamId, double>& p) const {
        ObsSamples out; 
        out.reserve(cfg_.draws);

        // auto start_smpl = std::chrono::steady_clock::now();
        NuisanceSamples samples = sampler_.sample(cfg_.draws);
        // auto stop_smpl  = std::chrono::steady_clock::now();
        // auto time_sampling_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop_smpl - start_smpl).count();
        // LOG_INFO("Sampling distribution took", time_sampling_ms, "ms.");

        // auto start_pred = std::chrono::steady_clock::now();
        std::size_t i = 1;
        for (const auto& s : samples) {
            LOG_INFO("Sample", i, "of", cfg_.draws);
            auto res = model_->predict_optimized(p, s);
            auto unzipped_res = flatten(res);
            std::map<BinnedObservableId, double> value = zip(unzipped_res.ids, unzipped_res.vals);
            out.emplace_back(std::move(value));
            i++;
        }
        // auto stop_pred  = std::chrono::steady_clock::now();
        // auto time_prediction_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop_pred - start_pred).count();
        // LOG_INFO("Predicting for sampled nuisances took", time_prediction_ms, "ms.");

        return MCRealization {out, samples};
    }

    MCResult summarize(const std::map<ParamId, double>& p) const {
        auto smpl = sample_predictions(p);

        std::ofstream fs;
        fs.open("obs_samples.csv");

        for (auto &&[oid, v] : smpl.sampled_obss[0])
            fs << oid.str() << ',';
        fs << '\n';

        for (auto &&ovec : smpl.sampled_obss) {
            for (auto &&[oid, v] : ovec)
                fs << v << ',';
            fs << '\n';
        }

        // auto start_sum = std::chrono::steady_clock::now();
        auto summary = gaussian_fit(smpl.sampled_obss, cfg_.skew_abs_threshold);
        // auto stop_sum  = std::chrono::steady_clock::now();
        // auto time_summarize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop_sum - start_sum).count();
        // LOG_INFO("Summarizing samples took", time_summarize_ms, "ms.");

        const auto covariance_ids = covariance_ids_from_first_sample(smpl.sampled_obss);
        auto covariance = covariance_from_obs_samples(
            smpl.sampled_obss,
            covariance_ids,
            cfg_.covariance_ridge_rel,
            cfg_.covariance_ridge_abs
        );

        return MCResult {smpl, summary, covariance};
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