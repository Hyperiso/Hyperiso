#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <array>

#include "Series.h"
#include "CSVOptions.h"

/**
 * @file DataFrame.h
 * @brief Minimal tabular container with column-wise storage and CSV I/O.
 *
 * This header provides a lightweight DataFrame type inspired by
 * common data-analysis libraries. It stores columns as Series<T>,
 * indexed by column name, and supports:
 *  - row/column access,
 *  - simple statistics and summary,
 *  - CSV export,
 *  - optional labeled index.
 */

/**
 * @brief Stream output helper for a vector of strings.
 *
 * Prints a vector as: [elem1, elem2, ...].
 */
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vec);

/**
 * @brief Stream output helper for a fixed-size array of 2 integers.
 *
 * Prints an array as: [a, b].
 */
std::ostream& operator<<(std::ostream& os, const std::array<int,2>& vec);

/**
 * @class DataFrame
 * @brief Simple column-oriented data structure with basic analysis utilities.
 *
 * The DataFrame class stores columns as Series<T>, where T is typically
 * int, double or std::string. Internally, columns are held in a map:
 *   - key   : column name (std::string),
 *   - value : std::shared_ptr<void> that actually points to Series<T>.
 *
 * Column types are tracked via a CSVOptions instance, which also
 * controls CSV export behavior (presence of index, etc.).
 *
 * Supported features include:
 *   - Adding columns and values dynamically,
 *   - Row-wise and label-based access (iat / at),
 *   - Simple statistical description for numeric columns,
 *   - Head / tail extraction,
 *   - CSV export via to_csv().
 */
class DataFrame {
private:
    /// Map from column name to a type-erased Series<T> (stored as shared_ptr<void>).
    std::unordered_map<std::string, std::shared_ptr<void>> columns_map;
    
    /// Current number of rows in the DataFrame.
    size_t nRows = 0;

    /// CSV-related options (column types, index presence, etc.).
    CSVOptions csvOptions;

    /// Internal size tracking (not exposed directly; shape is the public view).
    std::vector<int> size;
public:
    /// List of column names, in insertion order.
    std::vector<std::string> columns;

    /// Optional index labels for the rows (if set).
    std::vector<std::string> index;

    /**
     * @brief Shape of the DataFrame: {nRows, nCols}.
     *
     * shape[0] is the number of rows, shape[1] the number of columns.
     */
    std::array<int, 2> shape = {0,0};

    /**
     * @brief Default constructor.
     *
     * Constructs an empty DataFrame with an empty index.
     */
    DataFrame() : index(std::vector<std::string>()) {}

    /**
     * @brief Adds a new empty column of type T.
     *
     * Creates a new Series<T> with the given name and registers it in
     * the internal column map. The column type is also recorded in
     * csvOptions.columnTypes.
     *
     * @tparam T       Value type of the column.
     * @param colName  Name of the new column.
     */
    template < typename T>
    void addColumn(const std::string& colName);

    /**
     * @brief Appends a value to an existing column.
     *
     * The column must already exist and be of type T. If the column's size
     * exceeds the current number of rows, the DataFrame row count and shape
     * are updated accordingly.
     *
     * @tparam T       Value type of the column.
     * @param colName  Name of the target column.
     * @param value    Value to append.
     *
     * @throws std::invalid_argument if the column does not exist.
     */
    template <typename T>
    void addValueToColumn(const std::string& colName, const T& value);

    /**
     * @brief Positional access to a cell (row index, column name).
     *
     * @tparam T       Expected type of the column.
     * @param row      Zero-based row index.
     * @param colName  Name of the column.
     * @return Value at (row, colName).
     *
     * @throws std::invalid_argument if the column does not exist.
     * @throws std::out_of_range     if the row index is invalid.
     */
    template <typename T>
    T iat(size_t row, const std::string& colName) const;

    /**
     * @brief Label-based access to a cell (index label, column name).
     *
     * The DataFrame must have an index set (via setIndex). The method
     * searches the index for @p idxValue and returns the cell at the
     * corresponding row and column.
     *
     * @tparam T       Expected type of the column.
     * @param idxValue Row label to look up in the index.
     * @param colName  Name of the column.
     * @return Value at the matching row and column.
     *
     * @throws std::invalid_argument if the index is not set or if the label
     *         or column is not found.
     */
    template <typename T>
    T at(const std::string& idxValue, const std::string& colName) const;

    /**
     * @brief Accessor for a column as a modifiable Series<T>.
     *
     * If the column does not exist yet, it is created as an empty
     * Series<T> with the given name.
     *
     * @tparam T       Column value type.
     * @param colName  Name of the column.
     * @return Reference to the underlying Series<T>.
     */
    template <typename T>
    Series<T>& operator[](const std::string& colName);

    /**
     * @brief Prints the DataFrame to std::cout in a tabular format.
     *
     * Uses csvOptions.columnTypes to determine how to print each column
     * (int, double, std::string). Index labels are printed in the first
     * column if present.
     */
    void print() const;

    /**
     * @brief Returns the list of column names.
     */
    const std::vector<std::string>& getColumnNames() const;
    
    /**
     * @brief Returns the number of rows in the DataFrame.
     */
    size_t getRowCount() const;

    /**
     * @brief Returns a new DataFrame with the first n rows.
     *
     * Column metadata and csvOptions are propagated; index is truncated
     * consistently if present.
     *
     * @param n Number of rows to keep (clamped to available rows).
     * @return A new DataFrame containing the top n rows.
     */
    DataFrame head(size_t n = 5) const;

    /**
     * @brief Returns a new DataFrame with the last n rows.
     *
     * Column metadata and csvOptions are propagated; index is truncated
     * consistently if present.
     *
     * @param n Number of rows to keep (clamped to available rows).
     * @return A new DataFrame containing the bottom n rows.
     */
    DataFrame tail(size_t n = 5) const;

    /**
     * @brief Prints a statistical summary for a single numeric column.
     *
     * The summary includes:
     *   - count, mean, standard deviation,
     *   - min, quartiles (25%, 50%, 75%), and max.
     *
     * @tparam T       Numeric type of the column (e.g. int or double).
     * @param colName  Name of the column to describe.
     *
     * @throws std::invalid_argument if the column does not exist.
     * @throws std::runtime_error    if the series is empty.
     */
    template <typename T>
    void describeColumn(const std::string& colName) const;
    
    /**
     * @brief Prints a summary for all numeric columns.
     *
     * For each column whose type in csvOptions is int or double,
     * describeColumn() is called. Non-numeric columns are reported
     * as ignored.
     */
    void describe() const;

    /**
     * @brief Returns a copy of a column as a Series<T>.
     *
     * @tparam T       Value type of the column.
     * @param colName  Name of the column.
     * @return Copy of the Series<T> representing this column.
     *
     * @throws std::invalid_argument if the column does not exist or
     *         is not of type T.
     */
    template <typename T>
    Series<T> getColumn(const std::string& colName) const;

    /**
     * @brief Sets the index labels for all rows.
     *
     * The size of @p newIndex must match the current number of rows
     * (unless the DataFrame is still empty). The underlying Series<T>
     * of each column is updated to share the same index.
     *
     * @param newIndex Vector of labels to use as row index.
     *
     * @throws std::invalid_argument if the size does not match nRows.
     */

    void setIndex(const std::vector<std::string>& newIndex);

    /**
     * @brief Returns the current index labels.
     *
     * @return Reference to the index vector.
     *
     * @throws std::runtime_error if the index is not set.
     */
    const std::vector<std::string>& getIndex() const;

    /**
     * @brief Serializes the DataFrame to a CSV file.
     *
     * The header line contains the column names. If csvOptions.hasIndex
     * is true, an "Index" column is written first and each row begins
     * with its index (here: the integer row position).
     *
     * Column types come from csvOptions.columnTypes and must be known
     * (int, double, std::string) for all columns.
     *
     * @param filename Path to the target CSV file.
     *
     * @throws std::runtime_error on I/O errors or missing type info.
     */
    void to_csv(const std::string& filename);

    /**
     * @brief Internal helper to set CSV options from an external source.
     *
     * This copies the hasIndex flag and the columnTypes mapping, if
     * provided. Intended for use by CSV readers or DataFrame methods
     * that construct new DataFrames (e.g. head/tail).
     *
     * @param options Options to copy into this DataFrame.
     */
    void _set_csv_options(const CSVOptions& options);
};

#include "DataFrame.tpp"

#endif // DATAFRAME_H