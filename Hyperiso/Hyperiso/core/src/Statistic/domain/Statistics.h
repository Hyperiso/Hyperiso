#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>
#include <algorithm>
#include "Include.h"
#include "Math.h"

using Vec = std::vector<double>;
using Samples = std::vector<std::map<ObservableId, double>>; // shape: N x D (N samples of D-dim vector)

struct ColumnStats {
    double mean {0.0};
    double std_unbiased {0.0};
    double b1_skew {0.0};
    double std_p {0.0};
    double std_m {0.0};
    double mode {0.0};
};

inline std::map<ObservableId, ColumnStats> summarize_columns_obs(const Samples& S) {
    if (S.empty()) throw std::invalid_argument("No samples");

    const std::size_t N = S.size();
    const std::size_t Delta = std::floor(N * ERF_INV_RT2);
    const std::size_t D = S[0].size();

    for (const auto& v : S) {
        if (v.size()!=D) throw std::invalid_argument("Jagged samples");        
    }

    std::vector<ObservableId> ids;
    for (const auto& v : S[0]) {
        ids.push_back(v.first);
    }
    std::map<ObservableId, ColumnStats> out;

    std::vector<double> x;
    std::vector<double> w (N - Delta, 0.0);
    std::vector<double> theta (Delta, 0.0);
    for (size_t d = 0; d < D; d++) {
        for (const auto& v : S) x.push_back(v.at(ids[d]));
        std::sort(x.begin(), x.end());
        
        double mean {0.0};
        for (double x_i : x)  mean += x_i;
        mean /= static_cast<double>(N);
        
        double s = 0., m3 = 0.;
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

        for (std::size_t j = 0; j < N - Delta; j++)
            w[j] = x[j + Delta] - x[j];

        const std::size_t J = static_cast<std::size_t>(
            std::distance(w.begin(), std::min_element(w.begin(), w.end()))
        );

         const double width = w[J];
        if (width <= 0.0) {
            // degenerate, all equal in window
            const std::size_t K = J;
            out[ids[d]].mode  = x[K];
            out[ids[d]].std_m = 0.0;
            out[ids[d]].std_p = 0.0;
            continue;
        }

        for (std::size_t t = 0; t < Delta; t++) {
            const std::size_t k = J + t;
            theta[t] = std::abs(
                static_cast<double>(k) / static_cast<double>(N)
                - (x[k] - x[J]) / width
            );
        }

        const std::size_t tmin = static_cast<std::size_t>(
            std::distance(theta.begin(), std::min_element(theta.begin(), theta.end()))
        );
        const std::size_t K = J + tmin;

        out[ids[d]].mode  = x[K];
        out[ids[d]].std_m = x[K] - x[J];          // >= 0
        out[ids[d]].std_p = x[J + Delta] - x[K];  // >= 0

        x.clear();
    }
    
    return out;
}