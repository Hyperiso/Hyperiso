#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>
#include <algorithm>

#include "Include.h"
#include "Math.h"

/**
 * @file Statistics.h
 * @brief Statistical helpers for Monte Carlo observable samples.
 *
 * This header provides lightweight containers and inline utilities used to
 * summarize samples of binned observables and nuisance parameters.
 */

/** Convenient alias for a one-dimensional numeric vector. */
using Vec = std::vector<double>;

/** Observable samples stored as one observable-value map per Monte Carlo draw. */
using ObsSamples = std::vector<std::map<BinnedObservableId, double>>; // shape: N x D (N samples of D-dim vector)

/** Nuisance samples stored as one parameter-value map per Monte Carlo draw. */
using NuisanceSamples = std::vector<std::map<ParamId, double>>; // shape: N x D (N samples of D-dim vector)

/**
 * @struct ColumnStats
 * @brief Summary statistics for one sampled observable column.
 */
struct ColumnStats {
    double mean {0.0};          ///< Arithmetic mean of the sampled values.
    double std_unbiased {0.0};  ///< Unbiased standard deviation.
    double b1_skew {0.0};       ///< Moment skewness coefficient.
    double std_p {0.0};         ///< Right-side split standard deviation estimate.
    double std_m {0.0};         ///< Left-side split standard deviation estimate.
    double mode {0.0};          ///< Discrete mode-like central estimator.
};

/**
 * @brief Splits a vector around a reference value.
 *
 * Values lower than or equal to @p mu are returned in the first vector; values
 * above @p mu are returned in the second vector.
 *
 * @param x Input values.
 * @param mu Split threshold.
 * @return Pair @c {below_or_equal, above}.
 */
inline std::pair<std::vector<double>, std::vector<double>>
split_vector(const std::vector<double>& x, double mu)
{
    std::vector<double> below;
    std::vector<double> above;

    below.reserve(x.size());
    above.reserve(x.size());

    std::partition_copy(
        x.begin(), x.end(),
        std::back_inserter(below),
        std::back_inserter(above),
        [mu](double v) { return v <= mu; }
    );

    return {below, above};
}

/**
 * @brief Computes per-observable summary statistics from Monte Carlo samples.
 *
 * All sample maps must contain the same number of observables.  The observable
 * ordering is taken from the first sample.  For each observable, the function
 * computes the mean, unbiased standard deviation, skewness, a discrete
 * mode-like estimator, and split standard deviations around that central value.
 *
 * @param S Observable samples with shape conceptually equal to @c N x D.
 * @return Map from observable id to summary statistics.
 *
 * @throws std::invalid_argument if @p S is empty or if the sample maps are
 *         jagged.
 */
inline std::map<BinnedObservableId, ColumnStats> summarize_columns_obs(const ObsSamples& S) {
    if (S.empty()) throw std::invalid_argument("No samples");

    const std::size_t N = S.size();
    const std::size_t D = S[0].size();

    for (const auto& v : S) {
        if (v.size() != D) throw std::invalid_argument("Jagged samples");
    }

    std::vector<BinnedObservableId> ids;
    for (const auto& v : S[0]) {
        ids.push_back(v.first);
    }

    std::map<BinnedObservableId, ColumnStats> out;

    std::vector<double> x;
    x.reserve(N);

    for (size_t d = 0; d < D; d++) {
        x.clear(); 
        for (const auto& v : S) {
            x.push_back(v.at(ids[d]));
        }
        std::sort(x.begin(), x.end());

        double mean = 0.0;
        for (double x_i : x) mean += x_i;
        mean /= static_cast<double>(N);

        double s = 0.0, m3 = 0.0;
        for (double x_i : x) {
            const double r = x_i - mean;
            s += r * r;
            m3 += r * r * r;
        }

        out[ids[d]].mean = mean;
        out[ids[d]].std_unbiased = std::sqrt(s / static_cast<double>(N - 1));

        const double m2 = s / static_cast<double>(N);
        if (m2 > 0.0) {
            const double m3bar = m3 / static_cast<double>(N);
            out[ids[d]].b1_skew = m3bar / std::pow(m2, 1.5);
        } else {
            out[ids[d]].b1_skew = 0.0;
        }

        auto obj = [x, N] (std::size_t i) {
            double s_m = 0;
            double s_p = 0;
            for (size_t j = 0; j < i; j++) s_m += (x[j] - x[i]) * (x[j] - x[i]);
            for (size_t j = i; j < N; j++) s_p += (x[j] - x[i]) * (x[j] - x[i]);
            return std::cbrt(s_m) + std::cbrt(s_p);
        };

        double mu_hat = 0;
        double min_obj = std::numeric_limits<double>::infinity();

        for (std::size_t i = 0; i < N; ++i) {
            double obj_i = obj(i);

            if (obj_i < min_obj) {
                min_obj = obj_i;
                mu_hat = x[i];
            }
        }

        auto [x_m, x_p] = split_vector(x, mu_hat);
        double s_m = 0;
        double s_p = 0;
        for (double x_i : x_m) s_m += (x_i - mu_hat) * (x_i - mu_hat);
        for (double x_i : x_p) s_p += (x_i - mu_hat) * (x_i - mu_hat);

        out[ids[d]].mode  = mu_hat;
        out[ids[d]].std_m = std::sqrt(min_obj / N) * std::cbrt(s_m);
        out[ids[d]].std_p = std::sqrt(min_obj / N) * std::cbrt(s_p);
    }

    return out;
}

#endif