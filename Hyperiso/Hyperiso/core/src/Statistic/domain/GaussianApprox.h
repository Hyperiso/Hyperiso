#pragma once
#include <vector>
#include <stdexcept>
#include "Statistics.h"

struct GaussianSummary {
    BinnedObservableId id;
    double mu {};           // population mean
    double sigma {};        // population std
    double sigma_p {};      // population right std
    double sigma_m {};      // population left std
    double mode {};         // population mode
    double skew {};         // sample skewness
    bool symmetric {};      
};

inline std::ostream& operator<<(std::ostream& os, const GaussianSummary& gs) {
    os << ObservableMapper::str(gs.id.s) << " = ";
    if (gs.symmetric) {
        os << std::setprecision(4) << gs.mu << " +- " << gs.sigma << " (skew = " << std::setprecision(2) << gs.skew << ")";
    } else {
        os << std::setprecision(4) << gs.mode << " + " << gs.sigma_p << " - " << gs.sigma_m << " (skew = " << std::setprecision(2) << gs.skew << ")";
    }
    return os;
}

inline std::vector<GaussianSummary> gaussian_fit(const ObsSamples& S, double skew_abs_threshold=0.2) {
    const auto cols = summarize_columns_obs(S);
    std::vector<GaussianSummary> out;
    for (const auto& col : cols) {
        out.push_back(GaussianSummary{
            col.first, 
            col.second.mean, 
            col.second.std_unbiased,
            col.second.std_p,
            col.second.std_m,
            col.second.mode,
            col.second.b1_skew,
            std::fabs(col.second.b1_skew) < skew_abs_threshold
        });
    }

    return out;
}