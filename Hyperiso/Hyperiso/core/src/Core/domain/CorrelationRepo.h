#ifndef CORRELATIONREPO_H
#define CORRELATIONREPO_H

#include "Math.h"
#include "Include.h"

/**
 * @struct CorrelationMatrixPair
 * @brief Stores a pair of correlation matrices (statistical and systematic)
 *        for a given key type.
 *
 * The template parameter @p T is the type used as keys in the correlation
 * matrices (e.g. ParamId or ObservableId). Internally this relies on
 * SparseMatrix<T>, which behaves like an associative container:
 *
 *   - `contains(key)` tests presence of an entry.
 *   - `at(key)`       retrieves the stored double value.
 *   - iteration       yields pairs `(std::pair<T,T>, double)`.
 *
 * Symmetry:
 *   The helper `emplace` below automatically inserts entries both for
 *   (a,b) and (b,a), enforcing a symmetric correlation structure.
 *
 * @tparam T Key type used to index the correlation entries.
 */
template<typename T>
struct CorrelationMatrixPair {
    SparseMatrix<T> stat;   ///< Statistical correlation matrix.
    SparseMatrix<T> syst;   ///< Systematic correlation matrix.

    /**
     * @brief Retrieves the correlation values for a given pair.
     *
     * If a given (T,T) pair is not present in one of the matrices, the
     * corresponding correlation component is assumed to be 0.0.
     *
     * @param key A pair of keys (T, T) representing two parameters/observables.
     * @return A pair `(stat_value, syst_value)` for the requested key.
     */
    std::pair<double, double> at(const std::pair<T, T>& key) const {
        double statv = stat.contains(key) ? stat.at(key) : 0;
        double systv = syst.contains(key) ? syst.at(key) : 0;
        return std::make_pair(statv, systv);
    };

    /**
     * @brief Inserts a new entry into both the statistical and systematic
     *        matrices, enforcing symmetry.
     *
     * This function:
     *   - stores (key.first,  key.second) → stat_val / syst_val
     *   - stores (key.second, key.first) → stat_val / syst_val
     *
     * so that the resulting correlation matrices are symmetric by
     * construction.
     *
     * @param key       The pair of keys (T, T).
     * @param stat_val  Statistical correlation value.
     * @param syst_val  Systematic correlation value.
     */
    void emplace(std::pair<T, T>&& key, double stat_val, double syst_val) {
        stat.emplace(key, stat_val);
        syst.emplace(key, syst_val);
        stat.emplace(std::make_pair(key.second, key.first), stat_val);
        syst.emplace(std::make_pair(key.second, key.first), syst_val);
    };

    /**
     * @brief Stream output operator for a CorrelationMatrixPair.
     *
     * The output iterates over all entries in the statistical matrix and
     * prints them together with their corresponding systematic values:
     *
     *   (i, j): rho_stat + rho_syst
     *
     * The format is intended for debugging/logging rather than for
     * machine parsing.
     *
     * @tparam U Key type.
     * @param os  Output stream.
     * @param cmp The CorrelationMatrixPair to print.
     * @return The output stream.
     */
    template<typename U>
    friend std::ostream& operator<<(std::ostream&, const CorrelationMatrixPair<U>&);
};


/**
 * @class CorrelationRepository
 * @brief Manages correlations between parameters and between observables.
 *
 * This class acts as a central repository of correlation information.
 * It stores:
 *   - Parameter-level correlations: CorrelationMatrixPair<ParamId>
 *   - Observable-level correlations: CorrelationMatrixPair<ObservableId>
 *
 * Typical usage:
 *   - Call `set_correlation_matrix(...)` once to install initial matrices.
 *   - Optionally call `merge_correlation_matrix(...)` to overlay new
 *     entries (overwriting any existing ones with the same key).
 *   - Query correlations via `get_correlation(...)` and
 *     `get_combined_correlation(...)`.
 *
 * Note: The repository expects its internal shared_ptr members to be set
 * before any `get_correlation` calls are made.
 */
class CorrelationRepository {
public:
    /**
     * @brief Retrieves the correlation between two parameters.
     *
     * @param p1 First parameter ID.
     * @param p2 Second parameter ID.
     * @return A pair `(rho_stat, rho_syst)` with statistical and systematic
     *         correlation respectively.
     *
     * @pre `parameter_correlations` must be non-null and properly initialized.
     */
    std::pair<double, double> get_correlation(ParamId p1, ParamId p2) const;

    /**
     * @brief Retrieves the correlation between two observables given as enum.
     *
     * This overload converts the enum values to ObservableId using
     * ObservableMapper::to_id, then delegates to the (ObservableId,ObservableId)
     * overload.
     *
     * @param o1 First observable (enum).
     * @param o2 Second observable (enum).
     * @return A pair `(rho_stat, rho_syst)`.
     */
    std::pair<double, double> get_correlation(Observables o1, Observables o2) const;

    /**
     * @brief Retrieves the correlation between two observables given as IDs.
     *
     * @param o1 First observable ID.
     * @param o2 Second observable ID.
     * @return A pair `(rho_stat, rho_syst)`.
     *
     * @pre `observable_correlations` must be non-null and properly initialized.
     */
    std::pair<double, double> get_correlation(ObservableId o1, ObservableId o2) const;

    /**
     * @brief Computes the combined correlation (quadratic sum) between two
     *        entities of the same type.
     *
     * The combined correlation is defined as:
     *
     *   sqrt( rho_stat^2 + rho_syst^2 )
     *
     * which is simply `std::hypot(rho_stat, rho_syst)`.
     *
     * @tparam T Entity type (e.g. ParamId or Observables).
     * @param p1 First entity.
     * @param p2 Second entity.
     * @return The combined correlation.
     */
    template<typename T>
    double get_combined_correlation(T p1, T p2) const {
        auto corr = get_correlation(p1, p2);
        return std::hypot(corr.first, corr.second);
    }
    
    /**
     * @brief Sets (replaces) the correlation matrix for parameters.
     *
     * @param correlation_matrices Shared pointer to the parameter
     *        correlation matrices.
     */
    void set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrices);

    /**
     * @brief Sets (replaces) the correlation matrix for observables.
     *
     * @param correlation_matrix Shared pointer to the observable
     *        correlation matrices.
     */
    void set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ObservableId>> correlation_matrix);

    /**
     * @brief Merges a new correlation matrix into the existing parameter
     *        correlations.
     *
     * For each entry in `correlation_matrix->stat` and `correlation_matrix->syst`,
     * any existing entry in `parameter_correlations` with the same key is
     * overwritten via `insert_or_assign`.
     *
     * @param correlation_matrix Shared pointer to the new parameter
     *        correlation matrix.
     *
     * @pre `parameter_correlations` must be non-null.
     */
    void merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrix);

    /**
     * @brief Merges a new correlation matrix into the existing observable
     *        correlations.
     *
     * Same semantics as the parameter version, but applied to
     * `observable_correlations`.
     *
     * @param correlation_matrix Shared pointer to the new observable
     *        correlation matrix.
     *
     * @pre `observable_correlations` must be non-null.
     */
    void merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ObservableId>> correlation_matrix);

    /**
     * @brief Prints the current content of parameter and observable correlations.
     *
     * Uses the LOGGER macros to print:
     *   - all stored parameter correlations
     *   - all stored observable correlations
     *
     * in a diagnostic-friendly format (using the overloaded operator<< on
     * CorrelationMatrixPair).
     */
    void print_content() const;

private:
    std::shared_ptr<CorrelationMatrixPair<ParamId>> parameter_correlations;         ///< Parameter correlation matrices.
    std::shared_ptr<CorrelationMatrixPair<ObservableId>> observable_correlations;   ///< Observable correlation matrices.
};

#endif // CORRELATIONREPO_H
