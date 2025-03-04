#ifndef __CORRELATIONREPO_H__
#define __CORRELATIONREPO_H__

#include "Math.h"
#include "General.h"

template<typename T>
struct CorrelationMatrixPair {
    SparseMatrix<T> stat;
    SparseMatrix<T> syst;

    std::pair<double, double> at(const std::pair<T, T>& key) const {
        double statv = stat.contains(key) ? stat.at(key) : 0;
        double systv = syst.contains(key) ? syst.at(key) : 0;
        return std::make_pair(statv, systv);
    };

    void emplace(const std::pair<T, T>& key, double stat_val, double syst_val) const {
        stat.emplace(key, stat_val);
        stat.emplace(key.swap({key.second, key.first}), stat_val);
        syst.emplace(key, syst_val);
        syst.emplace(key.swap({key.second, key.first}), syst_val);
    };

    template<typename U>
    friend std::ostream& operator<<(std::ostream&, const CorrelationMatrixPair<U>&);
};

class CorrelationRepository {
public:
    std::pair<double, double> get_correlation(ParamId p1, ParamId p2);
    std::pair<double, double> get_correlation(Observables o1, Observables o2);

    template<typename T>
    double get_combined_correlation(T p1, T p2);
    
    void set_correlation_matrix(const CorrelationMatrixPair<ParamId>&& correlation_matrices);
    void set_correlation_matrix(const CorrelationMatrixPair<Observables>&& correlation_matrix);
    void merge_correlation_matrix(const CorrelationMatrixPair<ParamId>&& correlation_matrix);
    void merge_correlation_matrix(const CorrelationMatrixPair<Observables>&& correlation_matrix);

    void print_content();

private:
    static CorrelationMatrixPair<ParamId> parameter_correlations;
    static CorrelationMatrixPair<Observables> observable_correlations;
};

#endif // __CORRELATIONREPO_H__
