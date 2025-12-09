#ifndef SERIES_H
#define SERIES_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <memory>

/**
 * @file Series.h
 * @brief Simple one-dimensional data container with basic statistics.
 *
 * This header defines a templated Series type that behaves similarly to
 * a minimal 1D column:
 *   - it stores values of type T,
 *   - optionally associates a string index with each value,
 *   - exposes convenience methods for statistics (min, max, mean, stddev,
 *     quartiles).
 */

/**
 * @class Series
 * @brief Generic one-dimensional container of values with optional string index.
 *
 * The Series<T> class is intended as a lightweight helper for handling
 * small numerical (or generic) datasets. It supports:
 *   - positional access via iat(),
 *   - label-based access via at() when an index is provided,
 *   - descriptive statistics: min, max, mean, standard deviation,
 *     and quartiles.
 *
 * The index is shared via a std::shared_ptr so that multiple Series
 * can share the same labels if needed.
 *
 * @tparam T Value type (typically arithmetic, but not strictly required
 *           for all operations except statistics).
 */
template <typename T>
class Series {
private:
    std::string name;                                   ///< Optional column name.
    std::vector<T> data;                                ///< Stored data values.
    std::shared_ptr<std::vector<std::string>> index;    ///< Optional string index.
public:
    /**
     * @brief Default constructor.
     *
     * Constructs an empty series with an empty shared index.
     */
    Series() : index(std::make_shared<std::vector<std::string>>()) {}

    /**
     * @brief Constructs a Series with a given column name.
     *
     * @param colName Human-readable name of the series (e.g. "mass", "q2").
     */
    Series(const std::string& colName);

    /**
     * @brief Appends a value at the end of the series (without index).
     *
     * The associated index entry, if any, is not modified. This is mainly
     * useful for purely positional series.
     *
     * @param value Value to append.
     */
    void add(const T& value);

    /**
     * @brief Appends a value with an associated index label.
     *
     * Both the data vector and the shared index are extended.
     *
     * @param value Value to append.
     * @param idx   Index label for this value.
     */
    void add(const T& value, const std::string& idx) {
        data.push_back(value);
        index->push_back(idx);
    }

    /**
     * @brief Positional access to the i-th element.
     *
     * @param idx Zero-based position in the series.
     * @return Value at the requested position.
     *
     * @throws std::out_of_range if idx >= size().
     */
    T iat(size_t idx) const;

    /**
     * @brief Label-based access to a value.
     *
     * Searches the index for the given label and returns the corresponding
     * data value.
     *
     * @param idx Index label (string).
     * @return Value at the matching index position.
     *
     * @throws std::invalid_argument if the index label is not found.
     */
    T at(const std::string& idx) const;

    /**
     * @brief Returns the number of elements stored in the series.
     */
    size_t size() const;

    /**
     * @brief Returns the series name.
     *
     * @return Reference to the internal name string.
     */
    const std::string& getName() const;

    /**
     * @brief Prints all stored values to std::cout.
     *
     * This is mainly intended for debugging or quick inspection.
     */
    void print() const;

    /**
     * @brief Returns the minimum value in the series.
     *
     * @return Minimum element.
     *
     * @throws std::runtime_error if the series is empty.
     */
    T min() const;

    /**
     * @brief Returns the maximum value in the series.
     *
     * @return Maximum element.
     *
     * @throws std::runtime_error if the series is empty.
     */
    T max() const;

    /**
     * @brief Computes the arithmetic mean of the series.
     *
     * @return Mean value as a double.
     *
     * @throws std::runtime_error if the series is empty.
     */
    double mean() const;

    /**
     * @brief Computes the (population) standard deviation of the series.
     *
     * The variance is computed as:
     *   var = (1/N) * sum_i (x_i - mean)^2
     *
     * @return Standard deviation as a double.
     *
     * @throws std::runtime_error if the series is empty.
     */
    double stddev() const;

    /**
     * @brief Returns the first, second (median), and third quartiles.
     *
     * The method sorts a copy of the data and picks elements at
     * 25%, 50% and 75% of the size.
     *
     * @return Vector {Q1, Q2, Q3}.
     *
     * @throws std::runtime_error if the series is empty.
     */
    std::vector<T> quartiles() const;

    /**
     * @brief Accessor to the shared index vector.
     *
     * @return Shared pointer to the index storage.
     */
    const std::shared_ptr<std::vector<std::string>>& getIndex() const {
        return index;
    }

    /**
     * @brief Returns a copy of the underlying data as a std::vector.
     */
    std::vector<T> to_vec() {
        return data;
    }

    /**
     * @brief Returns a vector of stringified values.
     *
     * Each element is converted to std::string via std::to_string.
     * Intended for simple export/inspection use cases.
     */
    std::vector<std::string> to_string_vec() {
        std::vector<std::string> res{};
        for (auto _ : data) {
            res.push_back(std::to_string(_));
        }
        return res;
    }
};

#include "Series.tpp"

#endif // SERIES_H
