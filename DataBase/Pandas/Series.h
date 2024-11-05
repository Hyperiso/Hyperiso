#ifndef SERIES_H
#define SERIES_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

template <typename T>
class Series {
private:
    std::string name;
    std::vector<T> data;
    std::shared_ptr<std::vector<std::string>> index;
public:
    Series() : index(std::make_shared<std::vector<std::string>>()) {}
    Series(const std::string& colName);

    void add(const T& value);

    void add(const T& value, const std::string& idx) {
        data.push_back(value);
        index->push_back(idx);
    }

    T iat(size_t idx) const;

    T at(const std::string& idx) const {
        auto it = std::find(index->begin(), index->end(), idx);
        if (it != index->end()) {
            size_t pos = std::distance(index->begin(), it);
            return data.at(pos);
        }
        throw std::invalid_argument("Index not found");
    }

    size_t size() const;
    const std::string& getName() const;
    void print() const;

    T min() const;
    T max() const;
    double mean() const;
    double stddev() const;

    std::vector<T> quartiles() const;

    const std::shared_ptr<std::vector<std::string>>& getIndex() const {
        return index;
    }

    std::vector<T> to_vec() {
        return data;
    }

    std::vector<std::string> to_string_vec() {
        std::vector<std::string> res{};
        for (auto _ : data) {
            res.push_back(std::to_string(_));
        }
        return res;
    }
};

#include "Series.tpp"

#endif // SERIES_H
