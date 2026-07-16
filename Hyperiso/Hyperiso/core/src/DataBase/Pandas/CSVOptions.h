
#ifndef CSVOPTIONS_H
#define CSVOPTIONS_H

#include <string>
#include <vector>
#include <typeindex>
#include <unordered_map>

/**
 * @file CSVOptions.h
 * @brief Configuration holder for CSV import/export of a DataFrame.
 *
 * This struct captures how a DataFrame should be serialized to / deserialized
 * from CSV: whether an index is present, the type of each column, and the
 * type used for the index.
 */

/**
 * @struct CSVOptions
 * @brief Options describing how to interpret or write CSV data.
 *
 * The CSVOptions structure is used by DataFrame to:
 *  - remember whether the CSV file contains an explicit index column,
 *  - store the type of each named column,
 *  - optionally store the type of the index.
 *
 * Types are represented using std::type_index and are typically
 * int, double or std::string for this project.
 */
struct CSVOptions {
    /// True if the CSV representation contains a separate index column.
    bool hasIndex = false;

    /**
     * @brief Mapping from column name to its declared C++ type.
     *
     * The type_index values are usually typeid(int), typeid(double),
     * or typeid(std::string). Missing entries may indicate that the
     * column type has not yet been inferred or set.
     */
    std::unordered_map<std::string, std::type_index> columnTypes;

    /**
     * @brief Type used for the DataFrame index.
     *
     * Default is std::string, which matches the usual label-based index.
     */

    std::type_index indexType = typeid(std::string);

    /**
     * @brief Access or insert the type entry for a given column name.
     *
     * If the column is not present in @ref columnTypes, it is inserted
     * with a default placeholder type (typeid(void)) and a reference
     * to that entry is returned.
     *
     * This is mainly a convenience operator to allow syntax like:
     * @code
     * CSVOptions opt;
     * opt["pt"] = typeid(double);
     * @endcode
     *
     * @param colName Name of the column.
     * @return Reference to the associated std::type_index entry.
     */
    std::type_index& operator[](const std::string& colName) {
        auto result = columnTypes.emplace(colName, std::type_index(typeid(void)));
        return result.first->second;
    }
};

#endif