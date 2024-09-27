#include "DataFrame.h"

template <typename T>
void DataFrame::addColumn(const std::string& colName) {
    auto series = std::make_shared<Series<T>>(colName);
    columns_map[colName] = series;
    columns.push_back(colName);
}

template <typename T>
void DataFrame::addValueToColumn(const std::string& colName, const T& value) {
    auto col = std::static_pointer_cast<Series<T>>(columns_map[colName]);
    col->add(value);
    if (col->size() > nRows) {
        nRows = col->size();
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
void DataFrame::setIndex(const Series<T>& newIndex) {
    if (newIndex.size() != nRows && nRows > 0) {
        throw std::invalid_argument("Index size does not match number of rows");
    }
    index = std::make_shared<Series<T>>(newIndex);
}

template <typename T>
const Series<T>& DataFrame::getIndex() const {
    if (!index) {
        throw std::runtime_error("Index is not set");
    }
    return *std::static_pointer_cast<Series< T>>(index);
}