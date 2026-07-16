#include "CSVReader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <typeindex>
#include <memory>
#include <iostream>
#include <limits>
#include <cctype>

std::string trim_str(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    return (first == std::string::npos || last == std::string::npos) ? "" : str.substr(first, last - first + 1);
}

template <typename T>
T convertValue(const std::string& value) {
    std::string cleaned_value = trim_str(value);

    if constexpr (std::is_same<T, int>::value) {
        return std::stoi(cleaned_value);
    } else if constexpr (std::is_same<T, double>::value) {
        try {
            return std::stod(cleaned_value);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error converting '" << cleaned_value << "' to double" << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: value out of range for double: " << cleaned_value << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
    } else if constexpr (std::is_same<T, std::string>::value) {
        return cleaned_value;
    } else {
        throw std::invalid_argument("Unsupported type");
    }
}

template <typename T>
void addValueToDataFrame(DataFrame& df, const std::string& colName, const std::string& value) {
    df.addValueToColumn<T>(colName, convertValue<T>(value));
}

DataFrame CSVReader::read_csv(const std::string& filename, CSVOptions options) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "bite" << std::endl;
        throw std::runtime_error("Could not open file");
    }

    DataFrame df;
    std::string line;
    bool header = true;
    std::vector<std::string> headers;
    int index_id = 0;
    std::vector<std::string> index{};

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        if (header) {
            while (std::getline(ss, value, ',')) {
                headers.push_back(trim_str(value));

                const auto& colName = headers.back();

                auto it = options.columnTypes.find(colName);
                if (it != options.columnTypes.end()) {
                    const auto& colType = it->second;
                    if (colType == typeid(int)) {
                        df.addColumn<int>(colName);
                    } else if (colType == typeid(double)) {
                        df.addColumn<double>(colName);
                    } else if (colType == typeid(std::string)) {
                        df.addColumn<std::string>(colName);
                    } else {
                        throw std::invalid_argument("Unsupported column type");
                    }
                } else {
                    df.addColumn<double>(colName);
                    options.columnTypes.emplace(colName, std::type_index(typeid(double)));
                }
            }
            header = false;
        } else {
            size_t colIdx = 0;
            while (std::getline(ss, value, ',')) {
                if (!options.hasIndex && colIdx == 0) {
                    index.emplace_back(std::to_string(index_id++));
                }
                if (options.hasIndex && colIdx == 0) {
                    index.emplace_back(trim_str(value));
                } else {
                    const auto& colName = headers[colIdx];

                    auto it = options.columnTypes.find(colName);
                    if (it != options.columnTypes.end()) {
                        const auto& colType = it->second;
                        if (colType == typeid(int)) {
                            addValueToDataFrame<int>(df, colName, value);
                        } else if (colType == typeid(double)) {
                            addValueToDataFrame<double>(df, colName, value);
                        } else if (colType == typeid(std::string)) {
                            addValueToDataFrame<std::string>(df, colName, value);
                        } else {
                            throw std::invalid_argument("Unsupported column type");
                        }
                    } else {
                        addValueToDataFrame<double>(df, colName, value);
                    }
                }
                colIdx++;
            }
        }
    }

    df.setIndex(index);
    file.close();
    df._set_csv_options(options);

    return df;
}

