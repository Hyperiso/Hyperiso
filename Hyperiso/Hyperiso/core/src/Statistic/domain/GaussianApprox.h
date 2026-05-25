#ifndef GAUSSIANAPPROX_H
#define GAUSSIANAPPROX_H

#include <vector>
#include <stdexcept>

#include "Statistics.h"

/**
 * @file GaussianApprox.h
 * @brief Utilities for summarizing observable samples with Gaussian approximations.
 *
 * This header provides a lightweight summary structure and a helper that fits a
 * Gaussian or split-Gaussian approximation to each observable column, depending
 * on the measured skewness.
 *
 * @see Statistics.h
 */

/**
 * @struct GaussianSummary
 * @brief Summary statistics for one binned observable distribution.
 *
 * The fields store both symmetric and asymmetric uncertainty information. The
 * @ref symmetric flag indicates which representation should be preferred when
 * reporting the approximation.
 */
struct GaussianSummary {
    BinnedObservableId id;  ///< Binned observable identifier.
    double mu {};           // population mean
    double sigma {};        ///< Symmetric population standard deviation.
    double sigma_p {};      ///< Right-side standard deviation for an asymmetric approximation.
    double sigma_m {};      ///< Left-side standard deviation for an asymmetric approximation.
    double mode {};         ///< Estimated population mode.
    double skew {};         ///< Sample skewness estimator.
    bool symmetric {};      ///< True when the distribution is treated as sufficiently symmetric.
};

/**
 * @brief Streams a human-readable Gaussian summary.
 *
 * Symmetric summaries are printed as `mu +- sigma`; asymmetric summaries are
 * printed as `mode + sigma_p - sigma_m`. In both cases, the skewness is included
 * for diagnostics.
 *
 * @param os Output stream.
 * @param gs Summary to print.
 *
 * @return Reference to @p os.
 */
inline std::ostream& operator<<(std::ostream& os, const GaussianSummary& gs) {
    os << ObservableMapper::str(gs.id.s) << "[" << gs.id.p.first << "," <<gs.id.p.second << "] = ";
    if (gs.symmetric) {
        os << std::setprecision(4) << gs.mu << " +- " << gs.sigma << " (skew = " << std::setprecision(2) << gs.skew << ")";
    } else {
        os << std::setprecision(4) << gs.mode << " + " << gs.sigma_p << " - " << gs.sigma_m << " (skew = " << std::setprecision(2) << gs.skew << ")";
    }
    return os;
}

/**
 * @brief Builds Gaussian or split-Gaussian summaries for observable samples.
 *
 * Each observable column is summarized through @ref summarize_columns_obs. A
 * column is marked as symmetric when the absolute skewness is below
 * @p skew_abs_threshold; otherwise, the asymmetric mode and left/right standard
 * deviations are kept for reporting.
 *
 * @param S Observable samples to summarize.
 * @param skew_abs_threshold Maximum absolute skewness accepted for the symmetric
 *        Gaussian approximation.
 *
 * @return One @ref GaussianSummary per observable column.
 */
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

#endif