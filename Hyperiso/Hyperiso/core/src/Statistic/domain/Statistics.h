#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>
#include <algorithm>
#include "Include.h"

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
        for (const auto& v : S) x.emplace_back(v.at(ids[d]));   // Retrieve sample values
        std::sort(x.begin(), x.end());
        
        double mean {0.0};
        for (double x_i : x)  mean += x_i;
        mean /= N;
        
        double s, m3;
        for (double x_i : x) {
            double r = x_i - mean;
            s += r * r;
            m3 += r * r * r;
        }

        for (size_t j = 0; j < N - Delta; j++) {
            w[j] = x[j + Delta] - x[j];
        }
        std::size_t J = std::distance(w.begin(), std::min_element(w.begin(), w.end()));

        for (size_t k = J; k < J + Delta; k++) {
            theta[k] = std::abs((double) k / N - (x[k] - x[J]) / w[J]);
        }
        std::size_t K = std::distance(theta.begin(), std::min_element(theta.begin(), theta.end()));

        out[ids[d]].mean = mean;
        out[ids[d]].std_unbiased = std::sqrt(s / (N - 1));
        out[ids[d]].b1_skew = m3 / (N * std::pow(s, 3));
        out[ids[d]].mode = x[K];
        out[ids[d]].std_m = x[K] - x[J];
        out[ids[d]].std_p = x[J + Delta] - x[K];
    }
    
    return out;
}