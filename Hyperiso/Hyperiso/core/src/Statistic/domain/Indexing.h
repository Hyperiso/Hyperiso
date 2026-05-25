#ifndef INDEXING_H
#define INDEXING_H

#include <map>
#include <vector>
#include <stdexcept>

#include "Matrix.h"
#include "Include.h"
#include "ObservableValue.h"
#include "BinnedObservableId.h"

/**
 * @file Indexing.h
 * @brief Helpers for converting between indexed maps and dense containers.
 *
 * The utilities in this file are used to keep the association between logical
 * ids and numerical values while switching between STL maps, vectors, and matrix
 * representations.
 */

/**
 * @struct UnzipResult1D
 * @brief Result of converting a one-dimensional indexed map to dense vectors.
 *
 * @tparam T Type of the ids.
 * @tparam U Type of the associated values.
 */
template <typename T, typename U>
struct UnzipResult1D {
    std::vector<T> ids;     ///< Ids in map iteration order.
    std::vector<U> vals;    ///< Values associated with @ref ids.
};

/**
 * @struct UnzipResult2D
 * @brief Result of converting a square indexed matrix to dense vectors.
 *
 * @tparam T Type of the row and column ids.
 */
template <typename T>
struct UnzipResult2D {
    std::vector<T> ids;                     ///< Row and column ids in map iteration order.
    std::vector<std::vector<double>> vals;  ///< Dense square matrix values.
};

/**
 * @brief Builds a map from parallel id and value vectors.
 *
 * @tparam T Id type.
 * @tparam U Value type.
 *
 * @param ids Ordered ids.
 * @param vals Values associated with @p ids.
 *
 * @return Map associating each id with the value at the same position.
 *
 * @throws std::invalid_argument if the two vectors do not have the same size.
 */
template<typename T, typename U>
std::map<T, U> zip(const std::vector<T>& ids, const std::vector<U>& vals) {
    if (ids.size() != vals.size())
        throw std::invalid_argument("Index and value vector sizes don't match.");

    std::map<T, U> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        indexed.emplace(ids.at(i), vals.at(i));
    }

    return indexed;    
}

/**
 * @brief Builds a square nested map from ids and a dense matrix.
 *
 * @tparam T Row and column id type.
 *
 * @param ids Row and column ids.
 * @param vals Dense square matrix represented as nested vectors.
 *
 * @return Nested map such that `out[row_id][col_id]` stores the corresponding
 *         matrix entry.
 *
 * @throws std::invalid_argument if @p vals is empty, non-square with respect to
 *         @p ids, or has incompatible dimensions.
 */
template<typename T>
std::map<T, std::map<T, double>> zip(const std::vector<T>& ids, const std::vector<std::vector<double>>& vals) {
    if (vals.empty())
        throw std::invalid_argument("No values to zip.");

    if (ids.size() != vals.size() || ids.size() != vals.at(0).size())
        throw std::invalid_argument("Index and value sizes don't match or value matrix is not square.");

    std::map<T, std::map<T, double>> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        std::map<T, double> row;
        for (size_t j = 0; j < ids.size(); j++) {
            row.emplace(ids[j], vals[i][j]);
        }   
        indexed.emplace(ids[i], std::move(row));
    }

    return indexed;    
}

/**
 * @brief Flattens observable values into binned observable ids and values.
 *
 * Each input observable can contain several binned values. The output contains
 * one @ref BinnedObservableId per value, using `(0, 0)` as the bin range when an
 * observable value has no explicit bin.
 *
 * @param indexed Map from observable id to one or more observable values.
 *
 * @return Parallel vectors of binned ids and numerical values.
 */
inline UnzipResult1D<BinnedObservableId, double> flatten(std::map<ObservableId, std::vector<ObservableValue>> indexed) {
    std::vector<BinnedObservableId> ids;
    std::vector<double> vals;

    for (auto& [id, obs_values]: indexed) {
        for (auto& obs_val : obs_values) {
            BinnedObservableId binned_id;
            binned_id.s = id;
            binned_id.p = obs_val.bin.value_or(std::pair<double, double> {0., 0.});
            ids.emplace_back(binned_id);
            vals.emplace_back(obs_val.value);
        }
    }

    return UnzipResult1D {ids, vals};
}

/**
 * @brief Builds a square nested map from ids and a @ref RealMatrix.
 *
 * @tparam T Row and column id type.
 *
 * @param ids Row and column ids.
 * @param vals Dense square matrix.
 *
 * @return Nested map such that `out[row_id][col_id]` stores the corresponding
 *         matrix entry.
 *
 * @throws std::invalid_argument if @p vals is not square with respect to @p ids.
 */
template<typename T>
std::map<T, std::map<T, double>> zip(const std::vector<T>& ids, const RealMatrix& vals) {
    if (ids.size() != vals.rows() || ids.size() != vals.cols())
        throw std::invalid_argument("Index and value sizes don't match or value matrix is not square.");

    std::map<T, std::map<T, double>> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        std::map<T, double> row;
        for (size_t j = 0; j < ids.size(); j++) {
            row.emplace(ids[j], vals.unchecked_at(i, j));
        }   
        indexed.emplace(ids[i], std::move(row));
    }

    return indexed;  
}

/**
 * @brief Splits a map into parallel id and value vectors.
 *
 * @tparam T Id type.
 * @tparam U Value type.
 *
 * @param indexed Input map.
 *
 * @return Parallel vectors in map iteration order.
 */
template<typename T, typename U>
UnzipResult1D<T, U> unzip(const std::map<T, U>& indexed) {
    std::vector<T> ids;
    std::vector<U> vals;

    for (auto& [id, v]: indexed) {
        ids.emplace_back(id);
        vals.emplace_back(v);
    }

    return UnzipResult1D {ids, vals};
}

/**
 * @brief Splits a nested square map into ids and a dense matrix.
 *
 * @tparam T Row and column id type.
 *
 * @param indexed Nested square map.
 *
 * @return Id vector and dense matrix values in map iteration order.
 */
template<typename T>
UnzipResult2D<T> unzip(const std::map<T, std::map<T, double>>& indexed) {
    std::vector<T> ids;
    std::vector<std::vector<double>> vals;

    for (auto& [id, row]: indexed) {
        ids.emplace_back(id);
        std::vector<double> row_vals;
        for (auto& [_, v] : row) {
            row_vals.emplace_back(v);
        }
        vals.emplace_back(std::move(row_vals));
    }

    return UnzipResult2D {ids, vals};
}

#endif // INDEXING_H
