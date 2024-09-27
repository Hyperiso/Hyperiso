#include "Series.h"
#include <numeric>
#include <limits>
#include <stdexcept>
#include <algorithm>

template <typename T>
Series<T>::Series(const std::string& colName) : name(colName) {}

template <typename T>
void Series<T>::add(const T& value) {
    data.push_back(value);
}

template <typename T>
T Series<T>::iat(size_t idx) const {
    if (idx >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[idx];
}

template <typename T>
size_t Series<T>::size() const {
    return data.size();
}

template <typename T>
const std::string& Series<T>::getName() const {
    return name;
}

template <typename T>
void Series<T>::print() const {
    for (const auto& value : data) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}


template <typename T>
T Series<T>::min() const {
    if (data.empty()) throw std::runtime_error("Series is empty");
    return *std::min_element(data.begin(), data.end());
}

template <typename T>
T Series<T>::max() const {
    if (data.empty()) throw std::runtime_error("Series is empty");
    return *std::max_element(data.begin(), data.end());
}

template <typename T>
double Series<T>::mean() const {
    if (data.empty()) throw std::runtime_error("Series is empty");
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}


template <typename T>
double Series<T>::stddev() const {
    if (data.empty()) {
        throw std::runtime_error("Empty series");
    }

    double meanValue = mean();
    double sumOfSquares = 0.0;

    for (const auto& value : data) {
        sumOfSquares += (value - meanValue) * (value - meanValue);
    }

    return std::sqrt(sumOfSquares / data.size());
}

template <typename T>
std::vector<T> Series<T>::quartiles() const {
    if (data.empty()) {
        throw std::runtime_error("Empty series");
    }

    std::vector<T> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    auto q1 = sortedData[sortedData.size() * 0.25];
    auto q2 = sortedData[sortedData.size() * 0.50];  // Median
    auto q3 = sortedData[sortedData.size() * 0.75];

    return {q1, q2, q3};
}