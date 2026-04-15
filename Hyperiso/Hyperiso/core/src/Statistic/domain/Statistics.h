#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>
#include <algorithm>
#include "Include.h"
#include "Math.h"

using Vec = std::vector<double>;
using ObsSamples = std::vector<std::map<BinnedObservableId, double>>; // shape: N x D (N samples of D-dim vector)
using NuisanceSamples = std::vector<std::map<ParamId, double>>; // shape: N x D (N samples of D-dim vector)

struct ColumnStats {
    double mean {0.0};
    double std_unbiased {0.0};
    double b1_skew {0.0};
    double std_p {0.0};
    double std_m {0.0};
    double mode {0.0};
};

// inline std::map<BinnedObservableId, ColumnStats> summarize_columns_obs(const ObsSamples& S) {
//     if (S.empty()) throw std::invalid_argument("No samples");

//     const std::size_t N = S.size();
//     const std::size_t Delta = std::floor(N * ERF_INV_RT2);
//     const std::size_t D = S[0].size();

//     for (const auto& v : S) {
//         if (v.size()!=D) throw std::invalid_argument("Jagged samples");        
//     }

//     std::vector<BinnedObservableId> ids;
//     for (const auto& v : S[0]) {
//         ids.push_back(v.first);
//     }
//     std::map<BinnedObservableId, ColumnStats> out;

//     std::vector<double> x;
//     std::vector<double> w (N - Delta, 0.0);
//     std::vector<double> theta (Delta, 0.0);
//     for (size_t d = 0; d < D; d++) {
//         for (const auto& v : S) x.push_back(v.at(ids[d]));
//         std::sort(x.begin(), x.end());
        
//         double mean {0.0};
//         for (double x_i : x)  mean += x_i;
//         mean /= static_cast<double>(N);
        
//         double s = 0., m3 = 0.;
//         for (double x_i : x) {
//             const double r = x_i - mean;
//             s += r * r;
//             m3 += r * r * r;
//         }

//         out[ids[d]].mean = mean;
//         out[ids[d]].std_unbiased = std::sqrt(s / static_cast<double>(N - 1));

//         const double m2 = s / static_cast<double>(N);
//         if (m2 > 0.0) {
//             const double m3bar = m3 / static_cast<double>(N);
//             out[ids[d]].b1_skew = m3bar / std::pow(m2, 1.5);
//         } else {
//             out[ids[d]].b1_skew = 0.0;
//         }

//         for (std::size_t j = 0; j < N - Delta; j++)
//             w[j] = x[j + Delta] - x[j];

//         const std::size_t J = static_cast<std::size_t>(
//             std::distance(w.begin(), std::min_element(w.begin(), w.end()))
//         );

//          const double width = w[J];
//         if (width <= 0.0) {
//             // degenerate, all equal in window
//             const std::size_t K = J;
//             out[ids[d]].mode  = x[K];
//             out[ids[d]].std_m = 0.0;
//             out[ids[d]].std_p = 0.0;
//             continue;
//         }

//         for (std::size_t t = 0; t < Delta; t++) {
//             const std::size_t k = J + t;
//             theta[t] = std::abs(
//                 static_cast<double>(k) / static_cast<double>(N)
//                 - (x[k] - x[J]) / width
//             );
//         }

//         const std::size_t tmin = static_cast<std::size_t>(
//             std::distance(theta.begin(), std::min_element(theta.begin(), theta.end()))
//         );
//         const std::size_t K = J + tmin;

//         out[ids[d]].mode  = x[K];
//         out[ids[d]].std_m = x[K] - x[J];          // >= 0
//         out[ids[d]].std_p = x[J + Delta] - x[K];  // >= 0

//         x.clear();
//     }
    
//     return out;
// }


// Version until 26-04-13 with old deterministic fit of half-moments 
// inline std::map<BinnedObservableId, ColumnStats> summarize_columns_obs(const ObsSamples& S) {
//     if (S.empty()) throw std::invalid_argument("No samples");

//     const std::size_t N = S.size();
//     const std::size_t Delta = std::floor(N * ERF_INV_RT2);
//     const std::size_t D = S[0].size();

//     for (const auto& v : S) {
//         if (v.size() != D) throw std::invalid_argument("Jagged samples");
//     }

//     std::vector<BinnedObservableId> ids;
//     for (const auto& v : S[0]) {
//         ids.push_back(v.first);
//     }

//     std::map<BinnedObservableId, ColumnStats> out;

//     std::vector<double> x;
//     x.reserve(N);

//     std::vector<double> w(N - Delta, 0.0);
//     std::vector<double> theta(Delta, 0.0);

//     for (size_t d = 0; d < D; d++) {
//         x.clear(); 
//         for (const auto& v : S) {
//             x.push_back(v.at(ids[d]));
//         }
//         std::sort(x.begin(), x.end());

//         double mean = 0.0;
//         for (double x_i : x) mean += x_i;
//         mean /= static_cast<double>(N);

//         double s = 0.0, m3 = 0.0;
//         for (double x_i : x) {
//             const double r = x_i - mean;
//             s += r * r;
//             m3 += r * r * r;
//         }

//         out[ids[d]].mean = mean;
//         out[ids[d]].std_unbiased = std::sqrt(s / static_cast<double>(N - 1));

//         const double m2 = s / static_cast<double>(N);
//         if (m2 > 0.0) {
//             const double m3bar = m3 / static_cast<double>(N);
//             out[ids[d]].b1_skew = m3bar / std::pow(m2, 1.5);
//         } else {
//             out[ids[d]].b1_skew = 0.0;
//         }

//         for (std::size_t j = 0; j < N - Delta; j++) {
//             w[j] = x[j + Delta] - x[j];
//         }

//         const std::size_t J = static_cast<std::size_t>(
//             std::distance(w.begin(), std::min_element(w.begin(), w.end()))
//         );

//         const double width = w[J];
//         if (width <= 0.0) {
//             const std::size_t K = J;
//             out[ids[d]].mode  = x[K];
//             out[ids[d]].std_m = 0.0;
//             out[ids[d]].std_p = 0.0;
//             continue;
//         }

//         for (std::size_t t = 0; t < Delta; t++) {
//             const std::size_t k = J + t;
//             theta[t] = std::abs(
//                 static_cast<double>(k) / static_cast<double>(N)
//                 - (x[k] - x[J]) / width
//             );
//         }

//         const std::size_t tmin = static_cast<std::size_t>(
//             std::distance(theta.begin(), std::min_element(theta.begin(), theta.end()))
//         );
//         const std::size_t K = J + tmin;

//         out[ids[d]].mode  = x[K];
//         out[ids[d]].std_m = x[K] - x[J];
//         out[ids[d]].std_p = x[J + Delta] - x[K];
//     }

//     return out;
// }

// New version of 26-04-13 with proper fitting of split-normal moments
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