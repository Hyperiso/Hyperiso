#ifndef GRADIENT_HELPER_H
#define GRADIENT_HELPER_H

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Math.h"
#include "IProfileableLikelihood.h"

/**
 * @file GradientHelper.h
 * @brief Numerical derivative and Laplace profiling helpers for nuisance parameters.
 *
 * The helpers in this header support fast profiling of nuisance parameters in
 * likelihood scans. They provide finite-difference derivatives, Hessian
 * regularization, stationarity diagnostics and a hybrid Laplace/Newton
 * refinement scheme.
 *
 * @see IProfileableLikelihood
 * @see Profiler
 */

/**
 * @struct EtaDerivatives
 * @brief First-order derivatives of the likelihood with respect to selected nuisance parameters.
 *
 * The vector @ref g_eta contains finite-difference derivatives of the NLL with
 * respect to nuisance parameters. The matrix @ref J_eta contains the Jacobian
 * of model predictions with respect to the same nuisance directions.
 */
struct EtaDerivatives {
    std::vector<double> g_eta;  ///< NLL gradient restricted to the selected nuisance directions.
    RealMatrix J_eta;           ///< Observable Jacobian with rows as observables and columns as nuisance directions.
};

/**
 * @struct LaplaceProfileComputation
 * @brief Result of a Laplace nuisance-profile computation.
 */
struct LaplaceProfileComputation {
    double nll_hat = 1e300;         ///< Profiled or approximate profiled NLL value.
    std::vector<double> eta_hat;    ///< Estimated profiled nuisance vector.
    bool ok = false;                ///< True when the computation produced a finite, usable result.
};

/**
 * @struct LaplaceProfileOptions
 * @brief Numerical controls for the hybrid Laplace/Newton nuisance profiler.
 */
struct LaplaceProfileOptions {
    // If |dNLL/deta_i| * sigma_i is above this value after the analytic step,
    // the direction is considered locally ill-behaved and receives Newton refinement.
    double stationarity_threshold = 5e-2;           ///< Threshold on \f$|\partial\mathrm{NLL}/\partial\eta_i|\,\sigma_i\f$ above which a direction is refined.

    // Hard cap to avoid accidentally making a large expensive correction problem.
    std::size_t max_refined_eta = 4;                ///< Maximum number of nuisance directions corrected by Newton refinement.

    // Number of cheap outer correction cycles.
    std::size_t max_refinement_iters = 2;           ///< Maximum number of outer correction cycles.

    // Damping line-search parameters for the Newton correction.
    std::size_t max_line_search_halvings = 8;       ///< Maximum number of backtracking halvings for a Newton step.
    double max_newton_step_in_sigma = 1.0;          ///< Maximum Newton displacement measured in nuisance standard deviations.

    double hessian_eig_floor_rel = 1e-8;            ///< Relative eigenvalue floor used when regularizing Hessians.
    bool use_direct_nll_for_final_value = false;    ///< If true, use the direct NLL at the final profiled point as the reported value.

    bool debug_refinement = false;                  ///< If true, print stationarity and refinement diagnostics.
    std::size_t debug_top_eta = 8;                  ///< Maximum number of nuisance directions shown in debug output.
    std::string debug_label {};                     ///< Optional label appended to debug messages.
};

/**
 * @brief Computes a robust central finite-difference step.
 *
 * The step scales with the magnitude of @p x and is optionally limited by a
 * parameter-specific step hint.
 *
 * @param x Point at which the derivative is evaluated.
 * @param step_hint Optional scale hint, typically a parameter standard deviation.
 *
 * @return Positive finite-difference step.
 */
static double fd_step(double x, double step_hint) {
    const double abs_x = std::abs(x);
    const double abs_hint = std::abs(step_hint);

    double h = 1e-5 * std::max(1.0, abs_x);

    if (std::isfinite(abs_hint) && abs_hint > 0.0) {
        h = std::min(h, 1e-2 * abs_hint);
    }

    return std::max(h, 1e-8);
}

/**
 * @brief Computes a finite-difference step for a nuisance parameter.
 *
 * The returned step is reduced when parameter limits are present so that the
 * symmetric finite-difference stencil remains inside the allowed interval.
 *
 * @param def Parameter definition containing step hints and optional limits.
 * @param x Current nuisance value.
 *
 * @return Positive finite-difference step compatible with the parameter limits.
 *
 * @throws std::runtime_error if no finite positive step can be constructed.
 */
static double eta_fd_step_with_limits(
    const fit_app::ParameterDefinition& def,
    double x
) {
    double h = fd_step(x, def.step_hint);

    if (def.limits.has_value()) {
        const auto [lo, hi] = *def.limits;
        const double room_minus = x - lo;
        const double room_plus = hi - x;

        const double max_symmetric_h = 0.45 * std::min(room_minus, room_plus);
        if (max_symmetric_h > 0.0 && std::isfinite(max_symmetric_h)) {
            h = std::min(h, max_symmetric_h);
        }
    }

    if (!(h > 0.0) || !std::isfinite(h)) {
        throw std::runtime_error("Invalid finite-difference step for nuisance derivative");
    }

    return h;
}

/**
 * @brief Builds the ordered list of all nuisance indices.
 *
 * @param eta_dim Number of nuisance parameters.
 *
 * @return Vector containing indices from 0 to @p eta_dim - 1.
 */
static std::vector<std::size_t> all_eta_indices(std::size_t eta_dim) {
    std::vector<std::size_t> out(eta_dim);
    for (std::size_t i = 0; i < eta_dim; ++i) {
        out[i] = i;
    }
    return out;
}

/**
 * @brief Computes the complement of a nuisance-index subset.
 *
 * @param eta_dim Total number of nuisance parameters.
 * @param excluded Indices to remove from the full set.
 *
 * @return Sorted vector containing all nuisance indices not listed in @p excluded.
 */
static std::vector<std::size_t> eta_index_complement(
    std::size_t eta_dim,
    const std::vector<std::size_t>& excluded
) {
    std::set<std::size_t> excluded_set(excluded.begin(), excluded.end());

    std::vector<std::size_t> out;
    out.reserve(eta_dim);

    for (std::size_t i = 0; i < eta_dim; ++i) {
        if (!excluded_set.contains(i)) {
            out.push_back(i);
        }
    }

    return out;
}

/**
 * @brief Tests whether a nuisance-index list contains a given index.
 *
 * @param indices Candidate index list.
 * @param idx Index to look for.
 *
 * @return True if @p idx is present in @p indices.
 */
static bool eta_index_contains(
    const std::vector<std::size_t>& indices,
    std::size_t idx
) {
    return std::find(indices.begin(), indices.end(), idx) != indices.end();
}

/**
 * @brief Extracts a principal submatrix using an explicit index list.
 *
 * @param M Source square matrix.
 * @param idx Row and column indices retained in the output.
 *
 * @return Matrix containing \f$M_{ij}\f$ for all selected index pairs.
 */
static RealMatrix principal_submatrix_by_indices(
    const RealMatrix& M,
    const std::vector<std::size_t>& idx
) {
    RealMatrix out(idx.size(), idx.size());

    for (std::size_t i = 0; i < idx.size(); ++i) {
        for (std::size_t j = 0; j < idx.size(); ++j) {
            out.at(i, j) = M.at(idx[i], idx[j]);
        }
    }

    return out;
}

/**
 * @brief Computes nuisance derivatives for a selected subset of directions.
 *
 * For each nuisance direction, the function evaluates central finite
 * differences of both the model prediction and the NLL.
 *
 * @param like Profileable likelihood to differentiate.
 * @param p Fixed parameter-of-interest vector.
 * @param eta_base Base nuisance point.
 * @param eta_indices Nuisance directions to differentiate.
 *
 * @return NLL gradient and observable Jacobian restricted to @p eta_indices.
 *
 * @throws std::runtime_error if indices are invalid or model dimensions change.
 */
static EtaDerivatives compute_eta_derivatives_subset(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta_base,
    const std::vector<std::size_t>& eta_indices
) {
    const auto defs = like.get_param_defs();

    const std::size_t p_dim = like.p_dimension();
    const std::size_t eta_dim = like.eta_dimension();

    const std::vector<double> f0 = like.predict(p, eta_base);
    const std::size_t n_obs = f0.size();

    EtaDerivatives out;
    out.g_eta.assign(eta_indices.size(), 0.0);
    out.J_eta = RealMatrix(n_obs, eta_indices.size());

    for (std::size_t col = 0; col < eta_indices.size(); ++col) {
        const std::size_t a = eta_indices[col];
        if (a >= eta_dim) {
            throw std::runtime_error("Eta derivative index out of range");
        }

        const std::size_t theta_index = p_dim + a;
        const auto& def = defs[theta_index];

        const double h = eta_fd_step_with_limits(def, eta_base[a]);

        std::vector<double> eta_plus = eta_base;
        std::vector<double> eta_minus = eta_base;

        eta_plus[a] += h;
        eta_minus[a] -= h;

        const std::vector<double> f_plus = like.predict(p, eta_plus);
        const std::vector<double> f_minus = like.predict(p, eta_minus);

        if (f_plus.size() != n_obs || f_minus.size() != n_obs) {
            throw std::runtime_error("Model prediction size changed during finite difference");
        }

        for (std::size_t k = 0; k < n_obs; ++k) {
            out.J_eta.at(k, col) = (f_plus[k] - f_minus[k]) / (2.0 * h);
        }

        const double nll_plus = like.nll_from_split(p, eta_plus);
        const double nll_minus = like.nll_from_split(p, eta_minus);

        out.g_eta[col] = (nll_plus - nll_minus) / (2.0 * h);
    }

    return out;
}

/**
 * @brief Computes nuisance derivatives for all nuisance directions.
 *
 * @param like Profileable likelihood to differentiate.
 * @param p Fixed parameter-of-interest vector.
 * @param eta0 Base nuisance point.
 *
 * @return NLL gradient and observable Jacobian for the full nuisance vector.
 */
static EtaDerivatives compute_eta_derivatives(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta0
) {
    return compute_eta_derivatives_subset(
        like,
        p,
        eta0,
        all_eta_indices(like.eta_dimension())
    );
}

/**
 * @brief Computes the finite-difference NLL gradient for selected nuisance parameters.
 *
 * @param like Profileable likelihood to differentiate.
 * @param p Fixed parameter-of-interest vector.
 * @param eta Nuisance point at which the gradient is evaluated.
 * @param eta_indices Nuisance directions included in the output.
 *
 * @return Gradient vector ordered as @p eta_indices.
 *
 * @throws std::runtime_error if a nuisance index is out of range.
 */
static std::vector<double> eta_gradient_nll_subset(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta,
    const std::vector<std::size_t>& eta_indices
) {
    const auto defs = like.get_param_defs();
    const std::size_t p_dim = like.p_dimension();
    const std::size_t eta_dim = like.eta_dimension();

    std::vector<double> g(eta_indices.size(), 0.0);

    for (std::size_t col = 0; col < eta_indices.size(); ++col) {
        const std::size_t a = eta_indices[col];
        if (a >= eta_dim) {
            throw std::runtime_error("eta_gradient_nll_subset: eta index out of range");
        }

        const auto& def = defs[p_dim + a];
        const double h = eta_fd_step_with_limits(def, eta[a]);

        std::vector<double> eta_plus = eta;
        std::vector<double> eta_minus = eta;
        eta_plus[a] += h;
        eta_minus[a] -= h;

        const double fp = like.nll_from_split(p, eta_plus);
        const double fm = like.nll_from_split(p, eta_minus);

        g[col] = (fp - fm) / (2.0 * h);
    }

    return g;
}

/**
 * @brief Computes the finite-difference NLL gradient for all nuisance parameters.
 *
 * @param like Profileable likelihood to differentiate.
 * @param p Fixed parameter-of-interest vector.
 * @param eta Nuisance point at which the gradient is evaluated.
 *
 * @return Full nuisance-gradient vector.
 */
static std::vector<double> eta_gradient_nll(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta
) {
    return eta_gradient_nll_subset(
        like,
        p,
        eta,
        all_eta_indices(like.eta_dimension())
    );
}

/**
 * @brief Multiplies a matrix by a vector.
 *
 * @param M Matrix operand.
 * @param v Vector operand.
 *
 * @return Product \f$M v\f$.
 *
 * @throws std::runtime_error if dimensions are incompatible.
 */
static std::vector<double> matvec(const RealMatrix& M, const std::vector<double>& v) {
    if (M.cols() != v.size()) {
        throw std::runtime_error("matvec: dimension mismatch");
    }

    std::vector<double> out(M.rows(), 0.0);

    for (std::size_t i = 0; i < M.rows(); ++i) {
        for (std::size_t j = 0; j < M.cols(); ++j) {
            out[i] += M.at(i, j) * v[j];
        }
    }

    return out;
}

/**
 * @brief Computes the Euclidean dot product of two vectors.
 *
 * @param a First vector.
 * @param b Second vector.
 *
 * @return Dot product \f$a^\top b\f$.
 *
 * @throws std::runtime_error if vector sizes differ.
 */
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("dot: dimension mismatch");
    }

    double out = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i) {
        out += a[i] * b[i];
    }

    return out;
}

/**
 * @brief Symmetrizes a square matrix after validating all entries.
 *
 * The function averages each matrix entry with its transpose counterpart and
 * then enforces exact symmetry. Non-finite values are rejected explicitly.
 *
 * @param H Matrix to validate and symmetrize.
 * @param label Human-readable label used in exception messages.
 *
 * @return Symmetric matrix.
 *
 * @throws std::runtime_error if @p H is non-square or contains non-finite entries.
 */
static RealMatrix force_symmetric_checked(
    const RealMatrix& H,
    const std::string& label = "matrix"
) {
    if (H.rows() != H.cols()) {
        throw std::runtime_error(label + " is not square");
    }

    RealMatrix sym(H.rows(), H.cols());

    double max_asym = 0.0;
    double max_abs = 0.0;

    for (std::size_t i = 0; i < H.rows(); ++i) {
        for (std::size_t j = 0; j < H.cols(); ++j) {
            const double a = H.at(i, j);
            const double b = H.at(j, i);

            if (!std::isfinite(a) || !std::isfinite(b)) {
                std::ostringstream oss;
                oss << label
                    << " contains non-finite value at ("
                    << i << "," << j << ") or transpose entry";
                throw std::runtime_error(oss.str());
            }

            const double v = 0.5 * (a + b);
            sym.at(i, j) = v;

            max_asym = std::max(max_asym, std::abs(a - b));
            max_abs = std::max(max_abs, std::max(std::abs(a), std::abs(b)));
        }
    }

    // Currently kept for easy debugging if this helper needs to report asymmetry.
    (void)max_asym;
    (void)max_abs;

    // Re-force exact symmetry entry-by-entry, in case RealMatrix operators
    // or storage order leave tiny inconsistencies.
    for (std::size_t i = 0; i < H.rows(); ++i) {
        for (std::size_t j = i + 1; j < H.cols(); ++j) {
            const double v = 0.5 * (sym.at(i, j) + sym.at(j, i));
            sym.at(i, j) = v;
            sym.at(j, i) = v;
        }
    }

    return sym;
}

/**
 * @brief Regularizes a symmetric matrix into a numerically positive-definite matrix.
 *
 * Eigenvalues below a relative floor are lifted while preserving the
 * eigenvectors, then exact symmetry is restored.
 *
 * @param H Candidate Hessian matrix.
 * @param rel_floor Relative eigenvalue floor with respect to the largest
 *                  positive eigenvalue.
 * @param label Human-readable label used in exception messages.
 *
 * @return Symmetric positive-definite regularized matrix.
 *
 * @throws std::runtime_error if eigendecomposition fails or no positive
 *         eigenvalue is available.
 */
static RealMatrix regularize_spd_local(
    const RealMatrix& H,
    double rel_floor = 1e-10,
    const std::string& label = "Laplace Hessian"
) {
    RealMatrix sym = force_symmetric_checked(H, label);

    EigenSystem eig;
    try {
        eig = sym.eig();
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << label << " eig() failed after explicit symmetrization: "
            << e.what();
        throw std::runtime_error(oss.str());
    }

    double max_pos = 0.0;
    for (std::size_t i = 0; i < eig.D.rows(); ++i) {
        const double ev = eig.D.at(i, i);
        if (!std::isfinite(ev)) {
            throw std::runtime_error(label + " has non-finite eigenvalue");
        }
        if (ev > 0.0) {
            max_pos = std::max(max_pos, ev);
        }
    }

    if (!(max_pos > 0.0)) {
        throw std::runtime_error(label + " is not positive definite");
    }

    const double floor = std::max(1e-12, rel_floor * max_pos);

    RealMatrix Dreg(eig.D.rows(), eig.D.cols());

    for (std::size_t i = 0; i < eig.D.rows(); ++i) {
        Dreg.at(i, i) = std::max(eig.D.at(i, i), floor);
    }

    RealMatrix out = eig.P * Dreg * eig.P.transpose();

    // Force symmetry one final time.
    return force_symmetric_checked(out, label + " regularized");
}

/**
 * @brief Ranks nuisance directions by scaled stationarity violation.
 *
 * The ranking score is \f$|\partial\mathrm{NLL}/\partial\eta_i|\,\sigma_i\f$
 * when a finite scale is available, and the absolute gradient otherwise.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 * @param eta Nuisance point to diagnose.
 *
 * @return Pairs of score and nuisance index sorted from largest to smallest score.
 */
static std::vector<std::pair<double, std::size_t>> rank_eta_stationarity(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta
) {
    const std::size_t p_dim = like.p_dimension();
    const std::size_t eta_dim = like.eta_dimension();
    const auto defs = like.get_param_defs();

    const std::vector<double> g = eta_gradient_nll(like, p, eta);

    std::vector<std::pair<double, std::size_t>> ranked;
    ranked.reserve(eta_dim);

    for (std::size_t a = 0; a < eta_dim; ++a) {
        const double sigma = std::abs(defs[p_dim + a].step_hint);
        const double scaled = (sigma > 0.0 && std::isfinite(sigma))
            ? std::abs(g[a]) * sigma
            : std::abs(g[a]);

        if (std::isfinite(scaled)) {
            ranked.push_back({scaled, a});
        }
    }

    std::sort(
        ranked.begin(),
        ranked.end(),
        [](const auto& lhs, const auto& rhs) {
            return lhs.first > rhs.first;
        }
    );

    return ranked;
}

/**
 * @brief Selects nuisance directions requiring Newton refinement.
 *
 * Directions are ranked by scaled stationarity violation and filtered according
 * to @ref LaplaceProfileOptions::stationarity_threshold and
 * @ref LaplaceProfileOptions::max_refined_eta.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 * @param eta Current nuisance point.
 * @param options Refinement-selection options.
 *
 * @return Sorted nuisance-index list selected for refinement.
 */
static std::vector<std::size_t> select_nonstationary_eta_indices(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta,
    const LaplaceProfileOptions& options
) {
    std::vector<std::pair<double, std::size_t>> ranked =
        rank_eta_stationarity(like, p, eta);

    ranked.erase(
        std::remove_if(
            ranked.begin(),
            ranked.end(),
            [&](const auto& item) {
                return item.first <= options.stationarity_threshold;
            }
        ),
        ranked.end()
    );

    if (ranked.size() > options.max_refined_eta) {
        ranked.resize(options.max_refined_eta);
    }

    std::vector<std::size_t> out;
    out.reserve(ranked.size());

    for (const auto& [_, idx] : ranked) {
        out.push_back(idx);
    }

    std::sort(out.begin(), out.end());
    return out;
}

/**
 * @brief Formats a nuisance-index list for diagnostic output.
 *
 * @param indices Nuisance indices to format.
 *
 * @return Comma-separated list enclosed in braces.
 */
static std::string eta_index_list_string(
    const std::vector<std::size_t>& indices
) {
    std::ostringstream oss;
    oss << "{";
    for (std::size_t i = 0; i < indices.size(); ++i) {
        if (i) oss << ",";
        oss << indices[i];
    }
    oss << "}";
    return oss.str();
}

/**
 * @brief Prints a compact stationarity and refinement diagnostic summary.
 *
 * Output is emitted only when @ref LaplaceProfileOptions::debug_refinement is enabled.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 * @param eta Current nuisance point.
 * @param options Debug and refinement options.
 * @param iter Refinement iteration index.
 * @param direct_nll Direct NLL evaluated at @p eta.
 * @param bad Selected nuisance indices that will be refined.
 */
static void debug_print_stationarity_summary(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta,
    const LaplaceProfileOptions& options,
    std::size_t iter,
    double direct_nll,
    const std::vector<std::size_t>& bad
) {
    if (!options.debug_refinement) {
        return;
    }

    const auto defs = like.get_param_defs();
    const std::size_t p_dim = like.p_dimension();

    const auto ranked = rank_eta_stationarity(like, p, eta);
    const double max_scaled = ranked.empty() ? 0.0 : ranked.front().first;

    std::cout << "[LAPLACE WARMUP] " << options.debug_label
              << " iter=" << iter
              << " direct_nll=" << std::setprecision(12) << direct_nll
              << " max_scaled_grad=" << max_scaled
              << " threshold=" << options.stationarity_threshold
              << " status=" << (bad.empty() ? "PURE_LAPLACE_OK" : "REFINE")
              << " refine_eta=" << eta_index_list_string(bad)
              << std::endl;

    const std::size_t n_show = std::min(options.debug_top_eta, ranked.size());
    for (std::size_t k = 0; k < n_show; ++k) {
        const double scaled = ranked[k].first;
        const std::size_t a = ranked[k].second;
        const bool will_refine = eta_index_contains(bad, a);

        std::cout << "  [LAPLACE WARMUP] rank=" << k
                  << " eta_idx=" << a
                  << " theta_idx=" << (p_dim + a)
                  << " name=" << defs[p_dim + a].name
                  << " scaled_grad=" << std::setprecision(12) << scaled
                  << " action=" << (will_refine ? "REFINE" : "laplace-only")
                  << std::endl;
    }
}

/**
 * @brief Computes a finite-difference Hessian on selected nuisance directions.
 *
 * Diagonal terms are evaluated with a second-order central stencil and
 * off-diagonal terms with a mixed central stencil.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 * @param eta Nuisance point at which the Hessian is evaluated.
 * @param eta_indices Nuisance directions retained in the Hessian.
 *
 * @return Symmetric finite-difference Hessian restricted to @p eta_indices.
 */
static RealMatrix numerical_eta_hessian_subset(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta,
    const std::vector<std::size_t>& eta_indices
) {
    const std::size_t m = eta_indices.size();
    const auto defs = like.get_param_defs();
    const std::size_t p_dim = like.p_dimension();

    RealMatrix H(m, m);
    if (m == 0) {
        return H;
    }

    std::vector<double> h(m, 0.0);
    for (std::size_t c = 0; c < m; ++c) {
        const std::size_t a = eta_indices[c];
        h[c] = eta_fd_step_with_limits(defs[p_dim + a], eta[a]);
    }

    const double f0 = like.nll_from_split(p, eta);

    for (std::size_t ci = 0; ci < m; ++ci) {
        const std::size_t ai = eta_indices[ci];

        std::vector<double> ep = eta;
        std::vector<double> em = eta;
        ep[ai] += h[ci];
        em[ai] -= h[ci];

        const double fp = like.nll_from_split(p, ep);
        const double fm = like.nll_from_split(p, em);

        H.at(ci, ci) = (fp - 2.0 * f0 + fm) / (h[ci] * h[ci]);

        for (std::size_t cj = ci + 1; cj < m; ++cj) {
            const std::size_t aj = eta_indices[cj];

            std::vector<double> epp = eta;
            std::vector<double> epm = eta;
            std::vector<double> emp = eta;
            std::vector<double> emm = eta;

            epp[ai] += h[ci]; epp[aj] += h[cj];
            epm[ai] += h[ci]; epm[aj] -= h[cj];
            emp[ai] -= h[ci]; emp[aj] += h[cj];
            emm[ai] -= h[ci]; emm[aj] -= h[cj];

            const double fpp = like.nll_from_split(p, epp);
            const double fpm = like.nll_from_split(p, epm);
            const double fmp = like.nll_from_split(p, emp);
            const double fmm = like.nll_from_split(p, emm);

            const double hij = (fpp - fpm - fmp + fmm) / (4.0 * h[ci] * h[cj]);

            H.at(ci, cj) = hij;
            H.at(cj, ci) = hij;
        }
    }

    return force_symmetric_checked(H, "Newton finite-difference H_bad raw");
}

/**
 * @brief Clamps nuisance values to their configured parameter limits.
 *
 * Parameters without explicit limits are left unchanged.
 *
 * @param like Profileable likelihood providing parameter definitions.
 * @param eta Nuisance vector modified in place.
 */
static void clamp_eta_to_limits(
    const IProfileableLikelihood& like,
    std::vector<double>& eta
) {
    const auto defs = like.get_param_defs();
    const std::size_t p_dim = like.p_dimension();

    for (std::size_t a = 0; a < eta.size(); ++a) {
        const auto& def = defs[p_dim + a];
        if (def.limits.has_value()) {
            const auto [lo, hi] = *def.limits;
            eta[a] = std::clamp(eta[a], lo, hi);
        }
    }
}

/**
 * @brief Limits a Newton step in units of nuisance-parameter scale hints.
 *
 * The full correction vector is rescaled uniformly if any selected component
 * exceeds the configured maximum step in sigma units.
 *
 * @param like Profileable likelihood providing parameter definitions.
 * @param eta_indices Nuisance directions associated with @p delta.
 * @param delta Newton correction vector modified in place.
 * @param max_step_in_sigma Maximum allowed component-wise displacement in
 *                          sigma units.
 */
static void clamp_eta_step_in_sigmas(
    const IProfileableLikelihood& like,
    const std::vector<std::size_t>& eta_indices,
    std::vector<double>& delta,
    double max_step_in_sigma
) {
    if (!(max_step_in_sigma > 0.0) || !std::isfinite(max_step_in_sigma)) {
        return;
    }

    const auto defs = like.get_param_defs();
    const std::size_t p_dim = like.p_dimension();

    double scale = 1.0;

    for (std::size_t c = 0; c < eta_indices.size(); ++c) {
        const std::size_t a = eta_indices[c];
        const double sigma = std::abs(defs[p_dim + a].step_hint);

        if (sigma > 0.0 && std::isfinite(sigma)) {
            const double allowed = max_step_in_sigma * sigma;
            const double step = std::abs(delta[c]);
            if (step > allowed && step > 0.0) {
                scale = std::min(scale, allowed / step);
            }
        }
    }

    for (double& v : delta) {
        v *= scale;
    }
}

/**
 * @brief Applies the Laplace approximation to a subset of nuisance parameters.
 *
 * The method forms the approximate nuisance Hessian
 * \f$J_\eta^\top W_{\mathrm{obs}}J_\eta + W_\eta\f$, regularizes it, and
 * applies the standard quadratic-profile correction.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 * @param eta_base Expansion point for the nuisance vector.
 * @param laplace_eta_indices Nuisance directions treated by the Laplace step.
 * @param hessian_eig_floor_rel Relative eigenvalue floor for Hessian regularization.
 * @param hessian_label Label used in Hessian-related diagnostics.
 *
 * @return Approximate profiled NLL and profiled nuisance estimate.
 *
 * @throws std::runtime_error if the nuisance base point has the wrong dimension
 *         or required matrix operations fail.
 */
static LaplaceProfileComputation laplace_profile_eta_subset(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const std::vector<double>& eta_base,
    const std::vector<std::size_t>& laplace_eta_indices,
    double hessian_eig_floor_rel = 1e-10,
    const std::string& hessian_label = "Laplace H_eta"
) {
    if (eta_base.size() != like.eta_dimension()) {
        throw std::runtime_error("laplace_profile_eta_subset: eta_base has wrong dimension");
    }

    LaplaceProfileComputation out;
    out.eta_hat = eta_base;

    const double nll0 = like.nll_from_split(p, eta_base);

    if (laplace_eta_indices.empty()) {
        out.nll_hat = nll0;
        out.ok = std::isfinite(out.nll_hat);
        return out;
    }

    const std::vector<double> r0 = like.residuals(p, eta_base);

    const RealMatrix W_obs = force_symmetric_checked(
        like.observable_curvature(r0),
        "W_obs"
    );

    const RealMatrix W_eta_full = force_symmetric_checked(
        like.nuisance_curvature(eta_base),
        "W_eta full"
    );

    const RealMatrix W_eta = force_symmetric_checked(
        principal_submatrix_by_indices(W_eta_full, laplace_eta_indices),
        "W_eta subset"
    );

    const EtaDerivatives der = compute_eta_derivatives_subset(
        like,
        p,
        eta_base,
        laplace_eta_indices
    );

    const RealMatrix H_raw =
        der.J_eta.transpose() * W_obs * der.J_eta + W_eta;

    const RealMatrix H = regularize_spd_local(H_raw, hessian_eig_floor_rel, hessian_label);

    const RealMatrix H_inv = H.inv();

    const std::vector<double> Hinv_g = matvec(H_inv, der.g_eta);

    const double correction = 0.5 * dot(der.g_eta, Hinv_g);

    out.nll_hat = nll0 - correction;

    for (std::size_t col = 0; col < laplace_eta_indices.size(); ++col) {
        const std::size_t a = laplace_eta_indices[col];
        out.eta_hat[a] -= Hinv_g[col];
    }

    clamp_eta_to_limits(like, out.eta_hat);

    out.ok = std::isfinite(out.nll_hat);
    return out;
}

/**
 * @brief Profiles nuisance parameters with a hybrid Laplace/Newton procedure.
 *
 * The procedure first performs a full Laplace step, then detects directions
 * with large residual stationarity violations and refines them with damped
 * Newton corrections while retaining a Laplace approximation for the remaining
 * directions.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 * @param options Numerical and diagnostic options.
 *
 * @return Profiled NLL approximation and nuisance estimate.
 */
static LaplaceProfileComputation laplace_profile_eta_refined(
    const IProfileableLikelihood& like,
    const std::vector<double>& p,
    const LaplaceProfileOptions& options = {}
) {
    const std::size_t eta_dim = like.eta_dimension();

    LaplaceProfileComputation best =
        laplace_profile_eta_subset(
            like,
            p,
            like.central_eta(),
            all_eta_indices(eta_dim),
            options.hessian_eig_floor_rel,
            "Laplace initial H_eta"
        );

    if (!best.ok) {
        return best;
    }

    double best_direct = like.nll_from_split(p, best.eta_hat);

    for (std::size_t iter = 0; iter < options.max_refinement_iters; ++iter) {
        const std::vector<std::size_t> bad =
            select_nonstationary_eta_indices(like, p, best.eta_hat, options);

        debug_print_stationarity_summary(
            like,
            p,
            best.eta_hat,
            options,
            iter,
            best_direct,
            bad
        );

        if (bad.empty()) {
            break;
        }

        const std::vector<std::size_t> laplace_indices =
            eta_index_complement(eta_dim, bad);

        const std::vector<double> g_bad =
            eta_gradient_nll_subset(like, p, best.eta_hat, bad);

        RealMatrix H_bad =
            numerical_eta_hessian_subset(like, p, best.eta_hat, bad);

        try {
            H_bad = regularize_spd_local(
                H_bad,
                options.hessian_eig_floor_rel,
                "Newton refine H_bad"
            );
        } catch (const std::exception& e) {
            // The Newton correction block is locally non-convex or numerically unusable.
            // Do not invalidate the whole profile point; keep the best point found so far.
            if (options.debug_refinement) {
                std::cout << "[LAPLACE WARMUP] " << options.debug_label
                          << " iter=" << iter
                          << " newton_hessian=rejected"
                          << " reason=" << e.what()
                          << " refined_eta=" << eta_index_list_string(bad)
                          << std::endl;
            }
            break;
        }

        std::vector<double> delta = matvec(H_bad.inv(), g_bad);
        for (double& v : delta) {
            v = -v;
        }

        clamp_eta_step_in_sigmas(
            like,
            bad,
            delta,
            options.max_newton_step_in_sigma
        );

        bool accepted = false;
        LaplaceProfileComputation accepted_comp = best;
        double accepted_direct = best_direct;

        for (std::size_t ls = 0; ls <= options.max_line_search_halvings; ++ls) {
            const double alpha = std::ldexp(1.0, -static_cast<int>(ls));

            std::vector<double> eta_trial = best.eta_hat;
            for (std::size_t c = 0; c < bad.size(); ++c) {
                eta_trial[bad[c]] += alpha * delta[c];
            }

            clamp_eta_to_limits(like, eta_trial);

            LaplaceProfileComputation trial =
                laplace_profile_eta_subset(
                    like,
                    p,
                    eta_trial,
                    laplace_indices,
                    options.hessian_eig_floor_rel,
                    "Laplace refined subset H_eta"
                );

            if (!trial.ok) {
                continue;
            }

            const double direct_trial = like.nll_from_split(p, trial.eta_hat);

            if (std::isfinite(direct_trial) &&
                (!std::isfinite(best_direct) || direct_trial <= best_direct + 1e-10)) {
                accepted = true;
                accepted_comp = std::move(trial);
                accepted_direct = direct_trial;
                break;
            }
        }

        if (!accepted) {
            if (options.debug_refinement) {
                std::cout << "[LAPLACE WARMUP] " << options.debug_label
                          << " iter=" << iter
                          << " newton_step=rejected"
                          << " refined_eta=" << eta_index_list_string(bad)
                          << std::endl;
            }
            break;
        }

        if (options.debug_refinement) {
            std::cout << "[LAPLACE WARMUP] " << options.debug_label
                      << " iter=" << iter
                      << " newton_step=accepted"
                      << " direct_before=" << std::setprecision(12) << best_direct
                      << " direct_after=" << accepted_direct
                      << " refined_eta=" << eta_index_list_string(bad)
                      << std::endl;
        }

        best = std::move(accepted_comp);
        best_direct = accepted_direct;
    }

    if (options.use_direct_nll_for_final_value && std::isfinite(best_direct)) {
        best.nll_hat = best_direct;
    }

    best.ok = best.ok && std::isfinite(best.nll_hat);
    return best;
}

/**
 * @brief Convenience wrapper for nuisance profiling with default options.
 *
 * @param like Profileable likelihood.
 * @param p Fixed parameter-of-interest vector.
 *
 * @return Profiled NLL approximation and nuisance estimate.
 */
static LaplaceProfileComputation laplace_profile_eta(
    const IProfileableLikelihood& like,
    const std::vector<double>& p
) {
    return laplace_profile_eta_refined(like, p);
}

#endif
