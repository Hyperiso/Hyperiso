#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <fstream>
#include "Series.h"
#include "CSVOptions.h"
#include <array>

inline std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vec) {
    std::string result{};
    result.append("[");
    for (const auto& elem : vec) {
        result.append(elem + ", ");
    }
    result.replace(result.size()-2,2, "");
    result.append("]");
    os << result;
    return os;
}
inline std::ostream& operator<<(std::ostream& os, const std::array<int,2>& vec) {
    std::string result{};
    result.append("[");
    for (const auto& elem : vec) {
        result.append(std::to_string(elem) + ", ");
    }
    result.replace(result.size()-2,2, "");
    result.append("]");
    os << result;
    return os;
}
class DataFrame {
private:
    std::unordered_map<std::string, std::shared_ptr<void>> columns_map;
    
    size_t nRows = 0;

    CSVOptions csvOptions; 
    std::vector<int> size;
public:
    std::vector<std::string> columns;
    std::vector<std::string> index;
    std::array<int, 2> shape = {0,0};


    DataFrame() : index(std::vector<std::string>()) {}

    template < typename T>
    void addColumn(const std::string& colName);

    template <typename T>
    void addValueToColumn(const std::string& colName, const T& value);

    template <typename T>
    T iat(size_t row, const std::string& colName) const;

    template <typename T>
    T at(const std::string& idxValue, const std::string& colName) const;

    template <typename T>
    Series<T>& operator[](const std::string& colName);

    void print() const;

    const std::vector<std::string>& getColumnNames() const;
    
    size_t getRowCount() const;

    DataFrame head(size_t n = 5) const;
    DataFrame tail(size_t n = 5) const;

    template <typename T>
    void describeColumn(const std::string& colName) const;
    
    void describe() const;

    template <typename T>
    Series<T> getColumn(const std::string& colName) const;

    void setIndex(const std::vector<std::string>& newIndex);

    const std::vector<std::string>& getIndex() const;

    void to_csv(const std::string& filename);

    void _set_csv_options(CSVOptions options) {
        this->csvOptions = options;
    }
};

#include "DataFrame.tpp"

#endif // DATAFRAME_H