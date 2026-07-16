#include "DataFrame.h"
#include <iostream>

std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vec) {
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
std::ostream& operator<<(std::ostream& os, const std::array<int,2>& vec) {
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

void DataFrame::print() const {
    if (!index.empty()) {
        std::cout << "Index\t";
    }

    for (const auto& colName : columns) {
        std::cout << colName << "\t";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < nRows; ++i) {
        if (!index.empty()) {
            std::cout << index[i] << "\t";
        }

        for (const auto& colName : columns) {
            try {

                if (this->csvOptions.columnTypes.at(colName) == typeid(int)) {
                    std::cout << iat<int>(i, colName) << "\t";
                } else if (this->csvOptions.columnTypes.at(colName) == typeid(double)) {
                    std::cout << iat<double>(i, colName) << "\t";
                } else if (this->csvOptions.columnTypes.at(colName) == typeid(std::string)) {
                    std::cout << iat<std::string>(i, colName) << "\t";
                } else {
                    std::cout << "bug" << std::endl;
                    std::cout << "NaN\t";
                }
            } catch (const std::exception& e) {
                std::cout << e.what() << " " << colName << std::endl;
                std::cout << "ERROR" << std::endl;
                std::cout << "NaN\t";
            }
        }
        std::cout << std::endl;
    }
}


const std::vector<std::string>& DataFrame::getColumnNames() const {
    return columns;
}

size_t DataFrame::getRowCount() const {
    return nRows;
}

DataFrame DataFrame::head(size_t n) const {
    DataFrame df;
    n = std::min(n, nRows);

    for (const auto& colName : columns) {
        const std::type_index& colType = csvOptions.columnTypes.at(colName);

        if (colType == typeid(int)) {
            df.addColumn<int>(colName);
            for (size_t i = 0; i < n; ++i) {
                df.addValueToColumn<int>(colName, iat<int>(i, colName));
            }
        } else if (colType == typeid(double)) {
            df.addColumn<double>(colName);
            for (size_t i = 0; i < n; ++i) {
                df.addValueToColumn<double>(colName, iat<double>(i, colName));
            }
        } else if (colType == typeid(std::string)) {
            df.addColumn<std::string>(colName);
            for (size_t i = 0; i < n; ++i) {
                df.addValueToColumn<std::string>(colName, iat<std::string>(i, colName));
            }
        }
    }

    df._set_csv_options(this->csvOptions);

    if (!index.empty()) {
        std::vector<std::string> ind(index.begin(), index.begin() + n);
        df.setIndex(ind);
    }

    return df;
}

DataFrame DataFrame::tail(size_t n) const {
    DataFrame df;
    n = std::min(n, nRows);

    for (const auto& colName : columns) {
        const std::type_index& colType = csvOptions.columnTypes.at(colName);

        if (colType == typeid(int)) {
            df.addColumn<int>(colName);
            for (size_t i = nRows - n; i < nRows; ++i)
                df.addValueToColumn<int>(colName, iat<int>(i, colName));
        } else if (colType == typeid(double)) {
            df.addColumn<double>(colName);
            for (size_t i = nRows - n; i < nRows; ++i)
                df.addValueToColumn<double>(colName, iat<double>(i, colName));
        } else if (colType == typeid(std::string)) {
            df.addColumn<std::string>(colName);
            for (size_t i = nRows - n; i < nRows; ++i)
                df.addValueToColumn<std::string>(colName, iat<std::string>(i, colName));
        }
    }

    df._set_csv_options(this->csvOptions);

    if (!index.empty()) {
        std::vector<std::string> ind(index.end() - n, index.end());
        df.setIndex(ind);
    }

    return df;
}

void DataFrame::describe() const {
    std::cout << "Description of numerical columns :" << std::endl;
    
    for (const auto& colName : columns) {

        if (csvOptions.columnTypes.at(colName) == typeid(int)) {
            describeColumn<int>(colName);
        } else if (csvOptions.columnTypes.at(colName) == typeid(double)) {
            describeColumn<double>(colName);
        } else {
            std::cout << "Column '" << colName << "' is not numerical, then ignored." << std::endl;
        }
    }
}

void DataFrame::to_csv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Could not open file");

    const auto& cols = this->getColumnNames();
    size_t n = this->getRowCount();

    if (this->csvOptions.hasIndex) {
        file << "Index";
        for (const auto& c : cols) file << "," << c;
        file << "\n";
    } else {
        for (size_t j = 0; j < cols.size(); ++j) {
            file << cols[j];
            if (j + 1 < cols.size()) file << ",";
        }
        file << "\n";
    }

    for (size_t i = 0; i < n; ++i) {
        if (this->csvOptions.hasIndex) file << i << ",";

        for (size_t j = 0; j < cols.size(); ++j) {
            const auto& c = cols[j];

            auto it = this->csvOptions.columnTypes.find(c);
            if (it == this->csvOptions.columnTypes.end()) {
                throw std::runtime_error("Missing type info for column: " + c);
            }
            const std::type_index& ty = it->second;

            if (ty == typeid(int)) {
                file << this->iat<int>(i, c);
            } else if (ty == typeid(double)) {
                file << this->iat<double>(i, c);
            } else if (ty == typeid(std::string)) {
                file << this->iat<std::string>(i, c);
            } else {
                throw std::runtime_error("Unsupported column type for column: " + c);
            }

            if (j + 1 < cols.size()) file << ",";
        }
        file << "\n";
    }

    std::cout << "DataFrame successfully written to " << filename << std::endl;
}


const std::vector<std::string>& DataFrame::getIndex() const {
    
    if (index.empty()) {
        throw std::runtime_error("Index is not set");
    }
    return index;
}

void DataFrame::setIndex(const std::vector<std::string>& newIndex) {
    if (newIndex.size() != nRows && nRows > 0) {
        throw std::invalid_argument("Index size does not match number of rows");
    }
    index = newIndex;

    for (const auto& colName : columns) {
        auto it = csvOptions.columnTypes.find(colName);
        if (it == csvOptions.columnTypes.end()) continue;
        const std::type_index& t = it->second;

        if (t == typeid(int)) {
            auto s = std::static_pointer_cast<Series<int>>(columns_map.at(colName));
            *s->getIndex() = index;
        } else if (t == typeid(double)) {
            auto s = std::static_pointer_cast<Series<double>>(columns_map.at(colName));
            *s->getIndex() = index;
        } else if (t == typeid(std::string)) {
            auto s = std::static_pointer_cast<Series<std::string>>(columns_map.at(colName));
            *s->getIndex() = index;
        }
    }
}

void DataFrame::_set_csv_options(const CSVOptions& options) {
        this->csvOptions.hasIndex = options.hasIndex;

        if (!options.columnTypes.empty()) {
            this->csvOptions.columnTypes = options.columnTypes;
        }
    }