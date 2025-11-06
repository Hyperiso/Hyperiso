#pragma once
#include <vector>
#include <stdexcept>
#include "Statistics.h"


struct GaussianSummary { // per observable
    double mu{};
    double sigma{}; // unbiased std from MC
    double skew{}; // sample skewness
    bool approx_ok{};
};


inline std::vector<GaussianSummary> gaussian_or_warn(const Samples& S, double skew_abs_threshold=0.2) {
    const auto cols = summarize_columns(S);
    std::vector<GaussianSummary> out(cols.size());
    for (std::size_t i=0;i<cols.size();++i) {
        out[i] = GaussianSummary{cols[i].mean, cols[i].std_unbiased, cols[i].skewness,
                                std::fabs(cols[i].skewness) <= skew_abs_threshold};
    }
    return out;
}