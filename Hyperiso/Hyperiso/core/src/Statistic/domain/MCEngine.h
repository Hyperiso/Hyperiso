#pragma once
#include <vector>
#include <random>
#include "INuisanceSampler.h"
#include "ports/IModel.h"
#include "Statistics.h"
#include "GaussianApprox.h"
#include "Indexing.h"

struct MCConfig { std::size_t draws=10000; double skew_abs_threshold=0.2; };
struct MCRealization {
    ObsSamples sampled_obss;
    NuisanceSamples sampled_params;
};

struct MCResult {
    MCRealization mc_real;
    std::vector<GaussianSummary> summary;
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

        return MCResult {smpl, summary};
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