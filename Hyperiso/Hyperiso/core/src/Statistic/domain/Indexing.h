#ifndef __INDEXING_H__
#define __INDEXING_H__

#include <map>
#include <vector>
#include <stdexcept>
#include "Matrix.h"
#include "Include.h"
#include "ObservableValue.h"
#include "BinnedObservableId.h"

template <typename T, typename U>
struct UnzipResult1D {
    std::vector<T> ids;
    std::vector<U> vals;
};

template <typename T>
struct UnzipResult2D {
    std::vector<T> ids;
    std::vector<std::vector<double>> vals;
};

template<typename T, typename U>
std::map<T, U> zip(const std::vector<T>& ids, const std::vector<U>& vals) {
    if (ids.size() != vals.size())
        throw std::invalid_argument("Index and value vector sizes don't match.");

    std::map<T, U> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        indexed.emplace(ids.at(i), vals.at(i));
    }

    return indexed;    
}

template<typename T>
std::map<T, std::map<T, double>> zip(const std::vector<T>& ids, const std::vector<std::vector<double>>& vals) {
    if (vals.empty())
        throw std::invalid_argument("No values to zip.");

    if (ids.size() != vals.size() || ids.size() != vals.at(0).size())
        throw std::invalid_argument("Index and value sizes don't match or value matrix is not square.");

    std::map<T, std::map<T, double>> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        std::map<T, double> row;
        for (size_t j = 0; j < ids.size(); j++) {
            row.emplace(ids[j], vals[i][j]);
        }   
        indexed.emplace(ids[i], std::move(row));
    }

    return indexed;    
}

std::map<ObservableId, std::vector<ObservableValue>> zip(const std::vector<BinnedObservableId>& ids, const std::vector<double>& vals) {
    if (vals.empty())
        throw std::invalid_argument("No values to zip.");

    if (ids.size() != vals.size())
        throw std::invalid_argument("Index and value sizes don't match or value matrix is not square.");

    std::map<ObservableId, std::vector<ObservableValue>> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        auto bin = ids.at(i).p;
        if (indexed.contains(ids.at(i).s)) {
            if (fpeq(bin.first, bin.second)) {
                indexed.at(ids.at(i).s).emplace_back(ObservableValue(ids.at(i).s, vals.at(i)));
            } else {
                indexed.at(ids.at(i).s).emplace_back(ObservableValue(ids.at(i).s, vals.at(i), bin));
            }
        } else {
            if (fpeq(bin.first, bin.second)) {
                indexed.insert({ids.at(i).s, {ObservableValue(ids.at(i).s, vals.at(i))}});
            } else {
                indexed.insert({ids.at(i).s, {ObservableValue(ids.at(i).s, vals.at(i), bin)}});
            }
        }
    }

    return indexed;    
}

UnzipResult1D<BinnedObservableId, double> unzip(std::map<ObservableId, std::vector<ObservableValue>> indexed) {
    std::vector<BinnedObservableId> ids;
    std::vector<double> vals;

    for (auto& [id, obs_values]: indexed) {
        for (auto& obs_val : obs_values) {
            BinnedObservableId binned_id;
            binned_id.s = id;
            binned_id.p = obs_val.bin.value_or(std::pair<double, double> {0., 0.});
            ids.emplace_back(binned_id);
            vals.emplace_back(obs_val.value);
        }
    }

    return UnzipResult1D {ids, vals};
}

template<typename T>
std::map<T, std::map<T, double>> zip(const std::vector<T>& ids, const RealMatrix& vals) {
    if (ids.size() != vals.rows() || ids.size() != vals.cols())
        throw std::invalid_argument("Index and value sizes don't match or value matrix is not square.");

    std::map<T, std::map<T, double>> indexed;
    for (size_t i = 0; i < ids.size(); i++) {
        std::map<T, double> row;
        for (size_t j = 0; j < ids.size(); j++) {
            row.emplace(ids[j], vals.unchecked_at(i, j));
        }   
        indexed.emplace(ids[i], std::move(row));
    }

    return indexed;  
}

template<typename T, typename U>
UnzipResult1D<T, U> unzip(const std::map<T, U>& indexed) {
    std::vector<T> ids;
    std::vector<U> vals;

    for (auto& [id, v]: indexed) {
        ids.emplace_back(id);
        vals.emplace_back(v);
    }

    return UnzipResult1D {ids, vals};
}

template<typename T>
UnzipResult2D<T> unzip(const std::map<T, std::map<T, double>>& indexed) {
    std::vector<T> ids;
    std::vector<std::vector<double>> vals;

    for (auto& [id, row]: indexed) {
        ids.emplace_back(id);
        std::vector<double> row_vals;
        for (auto& [_, v] : row) {
            row_vals.emplace_back(v);
        }
        vals.emplace_back(std::move(row_vals));
    }

    return UnzipResult2D {ids, vals};
}

#endif // __INDEXING_H__
