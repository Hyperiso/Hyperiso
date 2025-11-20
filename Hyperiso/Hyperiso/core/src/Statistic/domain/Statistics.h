#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>
#include "Include.h"

using Vec = std::vector<double>;
using Samples = std::vector<std::map<ObservableId, double>>; // shape: N x D (N samples of D-dim vector)
using Samples_old = std::vector<Vec>;

struct ColumnStats {
    double mean{0.0};
    double std_unbiased{0.0};
    double skewness{0.0};
};


inline std::vector<ColumnStats> summarize_columns(const Samples_old& S) {
    if (S.empty()) throw std::invalid_argument("No samples");
    const std::size_t N = S.size();
    const std::size_t D = S[0].size();
    for (const auto& v : S) {

        if (v.size()!=D) throw std::invalid_argument("Jagged samples");
    }

    std::vector<ColumnStats> out(D);
    // means
    for (std::size_t d=0; d<D; ++d) {
        double m=0.0; for (const auto& v : S) m += v[d]; m/=static_cast<double>(N);
        out[d].mean = m;
    }
    // std + skew
    for (std::size_t d=0; d<D; ++d) {
        double m = out[d].mean;
        double s2=0.0, m3=0.0;
        for (const auto& v : S) {
            double c = v[d]-m; s2 += c*c; m3 += c*c*c;
        }
        // unbiased sample std (ddof=1)
        const double var_unb = s2 / static_cast<double>(N-1);
        out[d].std_unbiased = std::sqrt(var_unb);
        const double n = static_cast<double>(N);
        const double s = std::sqrt(s2/n); // population std
        if (s>0.0) {
            const double g1 = (n*n/((n-1)*(n-2))) * (m3/(n* s*s*s));
            out[d].skewness = g1;
        } else out[d].skewness = 0.0;
    }
    return out;
}

inline std::map<ObservableId, ColumnStats> summarize_columns_obs(const Samples& S) {
    if (S.empty()) throw std::invalid_argument("No samples");
    const std::size_t N = S.size();
    const std::size_t D = S[0].size();
    std::vector<ObservableId> ids(D);
    for (const auto& v : S[0]) {
        ids.push_back(v.first);
    }
    for (const auto& v : S) {
        if (v.size()!=D) throw std::invalid_argument("Jagged samples");
        
    }
    std::map<ObservableId, ColumnStats> out;
    // means
    for (std::size_t d=0; d<D; ++d) {
        double m=0.0; for (const auto& v : S) m += v.at(ids[d]); m/=static_cast<double>(N);
        out[ids[d]].mean = m;
    }
    // std + skew
    for (std::size_t d=0; d<D; ++d) {
        double m = out[ids[d]].mean;
        double s2=0.0, m3=0.0;
        for (const auto& v : S) {
            double c = v.at(ids[d])-m; s2 += c*c; m3 += c*c*c;
        }
        // unbiased sample std (ddof=1)
        const double var_unb = s2 / static_cast<double>(N-1);
        out[ids[d]].std_unbiased = std::sqrt(var_unb);
        const double n = static_cast<double>(N);
        const double s = std::sqrt(s2/n); // population std
        if (s>0.0) {
            const double g1 = (n*n/((n-1)*(n-2))) * (m3/(n* s*s*s));
            out[ids[d]].skewness = g1;
        } else out[ids[d]].skewness = 0.0;
    }
    return out;
}