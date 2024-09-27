#include "CSVReader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <typeindex>
#include <memory>
#include <iostream>

// Fonction pour convertir une valeur en fonction de son type
template <typename T>
T convertValue(const std::string& value) {
    if constexpr (std::is_same<T, int>::value) {
        return std::stoi(value);
    } else if constexpr (std::is_same<T, double>::value) {
        return std::stod(value);
    } else if constexpr (std::is_same<T, std::string>::value) {
        return value;
    } else {
        throw std::invalid_argument("Unsupported type");
    }
}

// Fonction générique pour ajouter une valeur dans une colonne
template <typename T>
void addValueToDataFrame(DataFrame& df, const std::string& colName, const std::string& value) {
    df.addValueToColumn<T>(colName, convertValue<T>(value));
}

DataFrame CSVReader::read_csv(const std::string& filename, const CSVOptions& options) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    DataFrame df;
    std::string line;
    bool header = true;
    std::vector<std::string> headers;

    // Index temporaire
    std::shared_ptr<void> index;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        if (header) {
            // Lecture des en-têtes de colonne
            while (std::getline(ss, value, ',')) {
                headers.push_back(value);
                if (!options.hasIndex || headers.size() > 1) {
                    // Ajouter une colonne seulement si ce n'est pas un index
                    const auto& colName = headers.back();
                    const auto& colType = options.columnTypes.at(colName);

                    if (colType == typeid(int)) {
                        df.addColumn<int>(colName);
                    } else if (colType == typeid(double)) {
                        df.addColumn<double>(colName);
                    } else if (colType == typeid(std::string)) {
                        df.addColumn<std::string>(colName);
                    } else {
                        throw std::invalid_argument("Unsupported column type");
                    }
                }
            }
            header = false;
        } else {
            size_t colIdx = 0;
            size_t indexPos = 0;
            std::vector<std::string> rowValues;

            // Parcours des valeurs de la ligne
            while (std::getline(ss, value, ',')) {
                if (options.hasIndex && colIdx == 0) {
                    // Gestion de l'index
                    if (options.indexType == typeid(int)) {
                        if (!index) index = std::make_shared<Series<int>>(Series<int>("Index"));
                        auto& idxSeries = *std::static_pointer_cast<Series<int>>(index);
                        idxSeries.add(convertValue<int>(value));
                    } else if (options.indexType == typeid(std::string)) {
                        if (!index) index = std::make_shared<Series<std::string>>(Series<std::string>("Index"));
                        auto& idxSeries = *std::static_pointer_cast<Series<std::string>>(index);
                        idxSeries.add(convertValue<std::string>(value));
                    } else {
                        throw std::invalid_argument("Unsupported index type");
                    }
                } else {
                    // Gestion des colonnes
                    const auto& colName = headers[colIdx];
                    const auto& colType = options.columnTypes.at(colName);

                    if (colType == typeid(int)) {
                        addValueToDataFrame<int>(df, colName, value);
                    } else if (colType == typeid(double)) {
                        addValueToDataFrame<double>(df, colName, value);
                    } else if (colType == typeid(std::string)) {
                        addValueToDataFrame<std::string>(df, colName, value);
                    } else {
                        throw std::invalid_argument("Unsupported column type");
                    }
                }
                colIdx++;
            }
        }
    }

    if (options.hasIndex) {
        // Assigner l'index si nécessaire
        if (options.indexType == typeid(int)) {
            df.setIndex(*std::static_pointer_cast<Series< int>>(index));
        } else if (options.indexType == typeid(std::string)) {
            df.setIndex(*std::static_pointer_cast<Series<std::string>>(index));
        }
    }

    file.close();
    df._set_csv_options(options);

    return df;
}
