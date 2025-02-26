#ifndef JSON_PARSER_H
#define JSON_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "Utils.h"

/**
 * @struct Value
 * @brief Represents a data entry with statistical and systematic errors.
 */
struct Value {
    std::string name;
    double central_value;
    double stat_error;
    double syst_error;
};

/**
 * @struct Correlation
 * @brief Represents a correlation between two data entries.
 */
struct Correlation {
    std::string name1;
    std::string name2;
    double value;
};


/**
 * @brief Removes leading and trailing whitespace from a string.
 * @param str The string to trim.
 * @return A new string with leading and trailing whitespace removed.
 */
std::string trim(const std::string& str);

/**
 * @brief Removes surrounding quotes from a string.
 * @param str The string from which quotes should be removed.
 * @return A new string without surrounding quotes.
 */
std::string remove_quotes(const std::string& str);

/**
 * @brief Reads a JSON file and extracts values and correlations.
 *
 * Parses the JSON file specified by `filename` to extract `Value` entries and
 * `Correlation` pairs, which are then stored in the provided vectors.
 *
 * @param filename The path to the JSON file to read.
 * @param values Vector to store extracted `Value` objects.
 * @param correlations Vector to store extracted `Correlation` objects.
 */
void read_json(const std::string& filename, std::vector<Value>& values, std::vector<Correlation>& correlations);

#endif // JSON_PARSER_H
