#include "MCEngine.h"
#include "StatisticProgress.h"
#include "ParameterRuntimeContext.h"

#include <algorithm>
#include <atomic>
#include <exception>
#include <fstream>
#include <mutex>
#include <thread>

RealMatrix symmetrize_covariance_matrix2(RealMatrix cov) {
    if (cov.rows() != cov.cols()) {
        throw std::runtime_error("symmetrize_covariance_matrix: covariance must be square");
    }

    for (std::size_t i = 0; i < cov.rows(); ++i) {
        for (std::size_t j = i + 1; j < cov.cols(); ++j) {
            const double a = cov.at(i, j);
            const double b = cov.at(j, i);

            if (!std::isfinite(a) || !std::isfinite(b)) {
                throw std::runtime_error("symmetrize_covariance_matrix: non-finite covariance entry");
            }

            const double v = 0.5 * (a + b);
            cov.at(i, j) = v;
            cov.at(j, i) = v;
        }

        if (!std::isfinite(cov.at(i, i))) {
            throw std::runtime_error("symmetrize_covariance_matrix: non-finite covariance diagonal");
        }
    }

    return cov;
}


RealMatrix inverse_covariance_with_ridge2(
    RealMatrix cov,
    double ridge_rel,
    double ridge_abs
) {
    cov = symmetrize_covariance_matrix2(std::move(cov));

    const std::size_t n = cov.rows();
    if (n != cov.cols()) {
        throw std::runtime_error("inverse_covariance_with_ridge: covariance must be square");
    }

    std::vector<double> sigma(n);

    for (std::size_t i = 0; i < n; ++i) {
        const double vii = cov.at(i, i);

        if (!std::isfinite(vii) || vii <= 0.0) {
            std::ostringstream oss;
            oss << "inverse_covariance_with_ridge: non-positive variance at i="
                << i << ", variance=" << vii;
            throw std::runtime_error(oss.str());
        }

        sigma[i] = std::sqrt(vii);
    }

    // Build dimensionless correlation matrix.
    RealMatrix corr(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            corr.at(i, j) = cov.at(i, j) / (sigma[i] * sigma[j]);
        }
    }

    corr = symmetrize_covariance_matrix2(std::move(corr));

    // Ridge in correlation space, dimensionless.
    const double ridge = std::max(ridge_rel, ridge_abs);

    for (std::size_t i = 0; i < n; ++i) {
        corr.at(i, i) += ridge;
    }

    corr = symmetrize_covariance_matrix2(std::move(corr));

    RealMatrix corr_inv = corr.inv();

    // Convert back: Cov^{-1} = D^{-1} Corr^{-1} D^{-1}
    RealMatrix cov_inv(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            cov_inv.at(i, j) = corr_inv.at(i, j) / (sigma[i] * sigma[j]);
        }
    }

    return symmetrize_covariance_matrix2(std::move(cov_inv));
}

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
    // for (std::size_t i = 0; i < D; ++i) {
    //     trace += std::max(0.0, cov.at(i, i));
    // }

    // const double scale = D > 0 ? trace / static_cast<double>(D) : 1.0;
    // const double ridge = std::max(ridge_abs, ridge_rel * std::max(scale, 1.0));

    // for (std::size_t i = 0; i < D; ++i) {
    //     cov.at(i, i) += ridge;
    // }

    // MCObservableCovariance out;
    // out.ids = ids;
    // out.mean = mean;
    // out.covariance = cov;
    // out.covariance_inv = cov.inv();
    // return out;
    for (std::size_t i = 0; i < D; ++i) {
        for (std::size_t j = i + 1; j < D; ++j) {
            const double v = 0.5 * (cov.at(i, j) + cov.at(j, i));
            cov.at(i, j) = v;
            cov.at(j, i) = v;
        }
    }

    MCObservableCovariance out;
    out.ids = ids;
    out.mean = mean;

    // Important: covariance brute, sans ridge physique.
    out.covariance = cov;

    // Option 1 : si covariance_inv n'est pas utilisée dans le fit MLE,
    // tu peux soit ne pas l'utiliser, soit la régulariser séparément.
    // out.covariance_inv = inverse_covariance_with_ridge2(
    //     cov,
    //     ridge_rel,
    //     ridge_abs
    // );
    out.covariance_inv = RealMatrix(D, D);
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

MCRealization MonteCarloEngine::sample_predictions_serial(const std::map<ParamId, double>& p) const {
    ObsSamples out;
    out.reserve(cfg_.draws);

    NuisanceSamples accepted_samples;
    accepted_samples.reserve(cfg_.draws);

    std::size_t accepted = 0;
    std::size_t failures = 0;
    std::size_t attempts = 0;
    StatisticProgressReporter progress(
        cfg_.print_progress,
        cfg_.draws,
        cfg_.progress_probe_draws,
        cfg_.progress_update_every
    );

    while (accepted < cfg_.draws) {
        ++attempts;

        std::map<ParamId, double> s = sampler_.sample();

        try {
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
            progress.accepted(accepted);

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

    progress.finish(accepted, attempts, failures);

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

MCRealization MonteCarloEngine::sample_predictions(const std::map<ParamId, double>& p) const {
    if (cfg_.n_threads <= 1 || cfg_.draws <= 1) {
        return sample_predictions_serial(p);
    }

    if (!model_->can_clone_for_worker()) {
        LOG_WARN(
            "Parallel MC requested with",
            cfg_.n_threads,
            "threads, but the model cannot be cloned for workers. Falling back to serial MC."
        );
        return sample_predictions_serial(p);
    }

    return sample_predictions_parallel(p);
}

MCRealization MonteCarloEngine::sample_predictions_parallel(const std::map<ParamId, double>& p) const {
    const std::size_t n_workers = std::max<std::size_t>(1, std::min(cfg_.n_threads, cfg_.draws));

    std::unique_ptr<IModelThreadGuard> decay_thread_guard;
    if (cfg_.force_decay_threads_to_one) {
        decay_thread_guard = model_->force_decay_threads(cfg_.forced_decay_threads);
    }

    struct WorkerOutput {
        ObsSamples obss;
        NuisanceSamples params;
        std::size_t failures = 0;
        std::size_t attempts = 0;
    };

    std::vector<WorkerOutput> worker_outputs(n_workers);
    std::vector<std::thread> workers;
    workers.reserve(n_workers);

    std::mutex sampler_mutex;
    std::mutex progress_mutex;
    std::mutex exception_mutex;
    std::exception_ptr first_exception = nullptr;

    std::atomic<std::size_t> accepted_total {0};
    std::atomic<std::size_t> failures_total {0};
    std::atomic<std::size_t> attempts_total {0};
    std::atomic<bool> stop {false};

    StatisticProgressReporter progress(
        cfg_.print_progress,
        cfg_.draws,
        cfg_.progress_probe_draws,
        cfg_.progress_update_every
    );

    auto set_exception_once = [&](std::exception_ptr eptr) {
        std::lock_guard<std::mutex> lock(exception_mutex);
        if (!first_exception) {
            first_exception = eptr;
        }
        stop.store(true, std::memory_order_release);
    };

    const std::size_t base_target = cfg_.draws / n_workers;
    const std::size_t remainder = cfg_.draws % n_workers;

    for (std::size_t worker_id = 0; worker_id < n_workers; ++worker_id) {
        const std::size_t target = base_target + (worker_id < remainder ? 1 : 0);

        workers.emplace_back([&, worker_id, target]() {
            auto& local = worker_outputs[worker_id];
            local.obss.reserve(target);
            local.params.reserve(target);

            try {
                ParameterRuntimeContext runtime_context;
                ScopedParameterRuntimeContext runtime_guard(runtime_context);

                auto worker_model = model_->clone_for_worker();
                if (!worker_model) {
                    throw std::runtime_error("MC worker could not clone model");
                }

                std::unique_ptr<IModelThreadGuard> worker_decay_thread_guard;
                if (cfg_.force_decay_threads_to_one) {
                    worker_decay_thread_guard = worker_model->force_decay_threads(cfg_.forced_decay_threads);
                }

                while (local.obss.size() < target && !stop.load(std::memory_order_acquire)) {
                    ++local.attempts;
                    ++attempts_total;

                    std::map<ParamId, double> s;
                    {
                        std::lock_guard<std::mutex> lock(sampler_mutex);
                        s = sampler_.sample();
                    }

                    try {
                        auto res = worker_model->predict_optimized(p, s);
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

                        local.obss.emplace_back(std::move(value));
                        local.params.emplace_back(std::move(s));

                        const std::size_t accepted_now =
                            accepted_total.fetch_add(1, std::memory_order_acq_rel) + 1;
                        {
                            std::lock_guard<std::mutex> lock(progress_mutex);
                            progress.accepted(accepted_now);
                        }
                    } catch (const std::exception& e) {
                        ++local.failures;
                        const std::size_t failures_now =
                            failures_total.fetch_add(1, std::memory_order_acq_rel) + 1;

                        LOG_WARN(
                            "Rejected MC nuisance sample",
                            failures_now,
                            "in worker",
                            worker_id,
                            "while trying to fill local accepted sample",
                            local.obss.size() + 1,
                            "of",
                            target,
                            ":",
                            e.what()
                        );

                        if (!cfg_.retry_failed_predictions ||
                            failures_now > cfg_.max_prediction_failures) {
                            set_exception_once(std::current_exception());
                            break;
                        }
                    } catch (...) {
                        ++local.failures;
                        const std::size_t failures_now =
                            failures_total.fetch_add(1, std::memory_order_acq_rel) + 1;

                        LOG_WARN(
                            "Rejected MC nuisance sample",
                            failures_now,
                            "with unknown exception in worker",
                            worker_id
                        );

                        if (!cfg_.retry_failed_predictions ||
                            failures_now > cfg_.max_prediction_failures) {
                            set_exception_once(std::current_exception());
                            break;
                        }
                    }
                }
            } catch (...) {
                set_exception_once(std::current_exception());
            }
        });
    }

    for (auto& worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }

    if (first_exception) {
        std::rethrow_exception(first_exception);
    }

    ObsSamples out;
    out.reserve(cfg_.draws);
    NuisanceSamples accepted_samples;
    accepted_samples.reserve(cfg_.draws);

    for (auto& local : worker_outputs) {
        for (auto& row : local.obss) {
            out.emplace_back(std::move(row));
        }
        for (auto& row : local.params) {
            accepted_samples.emplace_back(std::move(row));
        }
    }

    progress.finish(accepted_total.load(), attempts_total.load(), failures_total.load());

    if (failures_total.load() > 0) {
        LOG_WARN(
            "MC sampling finished with",
            failures_total.load(),
            "rejected nuisance samples over",
            attempts_total.load(),
            "attempts."
        );
    }

    return MCRealization{out, accepted_samples};
}

MCResult MonteCarloEngine::summarize(const std::map<ParamId, double>& p) const {
    auto smpl = sample_predictions(p);

    if (cfg_.write_samples_csv && !smpl.sampled_obss.empty()) {
        std::ofstream fs;
        fs.open(cfg_.samples_csv_path);

        for (auto &&[oid, v] : smpl.sampled_obss[0])
            fs << oid.str() << ',';
        fs << '\n';

        for (auto &&ovec : smpl.sampled_obss) {
            for (auto &&[oid, v] : ovec)
                fs << v << ',';
            fs << '\n';
        }
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