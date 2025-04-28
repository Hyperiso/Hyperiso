/**
 * @file CorrelationRepo.h
 * @brief Defines structures and classes to manage statistical and systematic correlations between parameters and observables.
 *
 * This file introduces:
 * - CorrelationMatrixPair: a structure to hold both statistical and systematic correlation matrices.
 * - CorrelationRepository: a class to access, set, and merge correlation data.
 */

#ifndef CORRELATIONREPO_H
#define CORRELATIONREPO_H

#include "Math.h"
#include "General.h"

/**
 * @struct CorrelationMatrixPair
 * @brief Stores a pair of correlation matrices (statistical and systematic) for a given type.
 *
 * @tparam T The type used as keys in the correlation matrices (e.g., ParamId or Observables).
 */
template<typename T>
struct CorrelationMatrixPair {
    SparseMatrix<T> stat;   ///< Statistical correlation matrix.
    SparseMatrix<T> syst;   ///< Systematic correlation matrix.

    /**
     * @brief Retrieves the correlation values for a given pair.
     *
     * @param key A pair of objects (T, T) representing two parameters/observables.
     * @return A pair containing (statistical value, systematic value).
     */
    std::pair<double, double> at(const std::pair<T, T>& key) const {
        double statv = stat.contains(key) ? stat.at(key) : 0;
        double systv = syst.contains(key) ? syst.at(key) : 0;
        return std::make_pair(statv, systv);
    };

    /**
     * @brief Inserts a new entry into both the statistical and systematic matrices.
     *
     * Automatically enforces symmetry (inserts both (A,B) and (B,A)).
     *
     * @param key The pair of keys (T, T).
     * @param stat_val Statistical correlation value.
     * @param syst_val Systematic correlation value.
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
     * Prints all correlations contained in the statistical matrix, along with systematic values.
     *
     * @tparam U Type of the key.
     * @param os Output stream.
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
 * Provides utilities to retrieve, set, merge, and print correlation matrices.
 */
class CorrelationRepository {
public:
    /**
     * @brief Retrieves the correlation between two parameters.
     *
     * @param p1 First parameter ID.
     * @param p2 Second parameter ID.
     * @return A pair containing (statistical correlation, systematic correlation).
     */
    std::pair<double, double> get_correlation(ParamId p1, ParamId p2) const;

    /**
     * @brief Retrieves the correlation between two observables.
     *
     * @param o1 First observable.
     * @param o2 Second observable.
     * @return A pair containing (statistical correlation, systematic correlation).
     */
    std::pair<double, double> get_correlation(Observables o1, Observables o2) const;

    /**
     * @brief Computes the combined correlation (hypotenuse) between two entities.
     *
     * @tparam T The type of entity (ParamId or Observables).
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
     * @brief Sets the correlation matrix for parameters.
     *
     * @param correlation_matrices Shared pointer to the parameter correlation matrices.
     */
    void set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrices);

    /**
     * @brief Sets the correlation matrix for observables.
     *
     * @param correlation_matrix Shared pointer to the observable correlation matrices.
     */
    void set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<Observables>> correlation_matrix);

    /**
     * @brief Merges a new correlation matrix into the existing parameter correlations.
     *
     * Existing entries will be overwritten if keys overlap.
     *
     * @param correlation_matrix Shared pointer to the new parameter correlation matrix.
     */
    void merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrix);

    /**
     * @brief Merges a new correlation matrix into the existing observable correlations.
     *
     * Existing entries will be overwritten if keys overlap.
     *
     * @param correlation_matrix Shared pointer to the new observable correlation matrix.
     */
    void merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<Observables>> correlation_matrix);

    /**
     * @brief Prints the current content of parameter and observable correlations.
     */
    void print_content() const;

private:
    std::shared_ptr<CorrelationMatrixPair<ParamId>> parameter_correlations;         ///< Parameter correlation matrices.
    std::shared_ptr<CorrelationMatrixPair<Observables>> observable_correlations;    ///< Observable correlation matrices.
};

#endif // CORRELATIONREPO_H
