#include "MCEngine.h"

MCObservableCovariance covariance_from_obs_samples(
    const ObsSamples& S,
    const std::vector<BinnedObservableId>& ids,
    double ridge_rel,
    double ridge_abs
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

std::vector<BinnedObservableId> covariance_ids_from_first_sample(
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

MCRealization MonteCarloEngine::sample_predictions(const std::map<ParamId, double>& p) const {
    ObsSamples out;
    out.reserve(cfg_.draws);

    NuisanceSamples accepted_samples;
    accepted_samples.reserve(cfg_.draws);

    std::size_t accepted = 0;
    std::size_t failures = 0;
    std::size_t attempts = 0;

    while (accepted < cfg_.draws) {
        ++attempts;

        std::map<ParamId, double> s = sampler_.sample();

        try {
            LOG_INFO("Sample", accepted + 1, "of", cfg_.draws);

            auto res = model_->predict_optimized(p, s);
            auto unzipped_res = flatten(res);
            std::map<BinnedObservableId, double> value =
                zip(unzipped_res.ids, unzipped_res.vals);

            bool finite = true;
            for (const auto& [oid, v] : value) {
                if (!std::isfinite(v)) {
                    finite = false;
                    break;
                }
            }

            if (!finite) {
                throw std::runtime_error("MC prediction contains non-finite observable");
            }

            out.emplace_back(std::move(value));
            accepted_samples.emplace_back(std::move(s));
            ++accepted;

        } catch (const std::exception& e) {
            ++failures;

            LOG_WARN(
                "Rejected MC nuisance sample",
                failures,
                "while trying to fill accepted sample",
                accepted + 1,
                "of",
                cfg_.draws,
                ":",
                e.what()
            );

            if (!cfg_.retry_failed_predictions ||
                failures > cfg_.max_prediction_failures) {
                throw;
            }

            continue;
        } catch (...) {
            ++failures;

            LOG_WARN(
                "Rejected MC nuisance sample",
                failures,
                "with unknown exception while trying to fill accepted sample",
                accepted + 1,
                "of",
                cfg_.draws
            );

            if (!cfg_.retry_failed_predictions ||
                failures > cfg_.max_prediction_failures) {
                throw;
            }

            continue;
        }
    }

    if (failures > 0) {
        LOG_WARN(
            "MC sampling finished with",
            failures,
            "rejected nuisance samples over",
            attempts,
            "attempts."
        );
    }

    return MCRealization{out, accepted_samples};
}

MCResult MonteCarloEngine::summarize(const std::map<ParamId, double>& p) const {
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

    auto summary = gaussian_fit(smpl.sampled_obss, cfg_.skew_abs_threshold);

    const auto covariance_ids = covariance_ids_from_first_sample(smpl.sampled_obss);
    auto covariance = covariance_from_obs_samples(
        smpl.sampled_obss,
        covariance_ids,
        cfg_.covariance_ridge_rel,
        cfg_.covariance_ridge_abs
    );

    return MCResult {smpl, summary, covariance};
}