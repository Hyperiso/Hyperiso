#ifndef CSVREADER_H
#define CSVREADER_H

#include <string>
#include <vector>
#include <typeindex>
#include <unordered_map>

#include "DataFrame.h"
#include "CSVOptions.h"

/**
 * @file CSVReader.h
 * @brief Helper for loading CSV files into a DataFrame.
 *
 * This header declares a small utility class that reads a CSV file and
 * constructs a @ref DataFrame according to a set of options
 * (column types, presence of an index, etc.).
 */

/**
 * @class CSVReader
 * @brief CSV → DataFrame loader with basic type handling.
 *
 * CSVReader provides a single high-level method @ref read_csv which:
 *   - reads a CSV file from disk,
 *   - uses the first line as a header (column names),
 *   - optionally interprets the first column as an index (see CSVOptions::hasIndex),
 *   - converts each cell to the column's declared type (int, double, std::string),
 *   - returns a fully constructed @ref DataFrame with correct shape, columns
 *     and index set.
 *
 * If a column type is not specified in CSVOptions::columnTypes, it defaults
 * to double and the option map is updated accordingly.
 */
class CSVReader {
public:
    /**
     * @brief Read a CSV file and construct a DataFrame.
     *
     * The file is expected to have a header line with column names.
     * Types and index behavior are controlled via the @p options argument:
     *
     *  - options.hasIndex:
     *      * true  → first column in the CSV is treated as an index column,
     *                values are stored in DataFrame::index and not as a data column.
     *      * false → rows are given a synthetic integer index ("0", "1", ...).
     *
     *  - options.columnTypes:
     *      * if a column name appears in this map, its values are converted
     *        to the corresponding type (int, double, std::string),
     *      * otherwise the column is assumed to be double.
     *
     * The returned DataFrame has its internal CSVOptions synchronized with
     * the effective options used during parsing.
     *
     * @param filename Path to the CSV file on disk.
     * @param options  Parsing options (column types, index presence).
     * @return A DataFrame containing the parsed CSV data.
     *
     * @throws std::runtime_error    If the file cannot be opened.
     * @throws std::invalid_argument If an unsupported column type is requested.
     */
    DataFrame read_csv(const std::string& filename, CSVOptions options = CSVOptions());
};

#endif // CSVREADER_H
