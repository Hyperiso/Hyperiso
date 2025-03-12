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

    void emplace(std::pair<T, T>&& key, double stat_val, double syst_val) {
        stat.emplace(key, stat_val);
        syst.emplace(key, syst_val);
        stat.emplace(std::make_pair(key.second, key.first), stat_val);
        syst.emplace(std::make_pair(key.second, key.first), syst_val);
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
    
    void set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrices);
    void set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<Observables>> correlation_matrix);
    void merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrix);
    void merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<Observables>> correlation_matrix);

    void print_content();

private:
    std::shared_ptr<CorrelationMatrixPair<ParamId>> parameter_correlations;
    std::shared_ptr<CorrelationMatrixPair<Observables>> observable_correlations;
};

#endif // __CORRELATIONREPO_H__
