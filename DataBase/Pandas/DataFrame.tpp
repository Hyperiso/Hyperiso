#include "DataFrame.h"

template <typename T>
void DataFrame::addColumn(const std::string& colName) {
    auto series = std::make_shared<Series<T>>(colName);
    columns_map[colName] = series;
    columns.push_back(colName);
    shape[1] +=1;
    csvOptions.columnTypes.emplace(colName, typeid(T));
}

template <typename T>
void DataFrame::addValueToColumn(const std::string& colName, const T& value) {
    auto col = std::static_pointer_cast<Series<T>>(columns_map[colName]);
    col->add(value);
    if (col->size() > nRows) {
        nRows = col->size();
        shape[0] = nRows;
    }
    if (csvOptions.columnTypes.find(colName) == csvOptions.columnTypes.end()) {
        csvOptions.columnTypes.emplace(colName, typeid(T));
    }
}

template <typename T>
T DataFrame::iat(size_t row, const std::string& colName) const {
    if (columns_map.find(colName) == columns_map.end()) {
        throw std::invalid_argument("Column not found");
    }
    auto col = std::static_pointer_cast<Series<T>>(columns_map.at(colName));
    return col->iat(row);
}

template <typename T>
T DataFrame::at(const std::string& idxValue, const std::string& colName) const {
    if (columns_map.find(colName) == columns_map.end()) {
        throw std::invalid_argument("Column not found");
    }
    auto col = std::static_pointer_cast<Series<T>>(columns_map.at(colName));
    
    return col->at(idxValue);
}

template <typename T>
Series<T> DataFrame::getColumn(const std::string& colName) const {
    if (columns_map.find(colName) == columns_map.end()) {
        throw std::invalid_argument("Column not found");
    }
    return *std::static_pointer_cast<Series<T>>(columns_map.at(colName));
}

template <typename T>
Series<T>& DataFrame::operator[](const std::string& colName) {
    if (columns_map.find(colName) == columns_map.end()) {
        addColumn<T>(colName);
    }
    return *std::static_pointer_cast<Series<T>>(columns_map[colName]);
}

template <typename T>
void DataFrame::describeColumn(const std::string& colName) const {
    // auto series = std::make_shared<Series<T>>(colName);
    const auto& series = getColumn<T>(colName);  // Récupérer la colonne

    std::cout << "Colonne: " << colName << std::endl;
    std::cout << "Count: " << series.size() << std::endl;
    std::cout << "Mean: " << series.mean() << std::endl;
    std::cout << "Std: " << series.stddev() << std::endl;
    std::cout << "Min: " << series.min() << std::endl;

    auto quartiles = series.quartiles();
    std::cout << "25%: " << quartiles[0] << std::endl;
    std::cout << "50%: " << quartiles[1] << std::endl;  // Median
    std::cout << "75%: " << quartiles[2] << std::endl;
    std::cout << "Max: " << series.max() << std::endl;
    std::cout << std::endl;
}

