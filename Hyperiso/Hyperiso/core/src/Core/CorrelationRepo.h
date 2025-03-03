#ifndef __CORRELATIONREPO_H__
#define __CORRELATIONREPO_H__

#include <variant>
#include "Math.h"
#include "General.h"

class CorrelationRepository {
public:
    enum class CorrelationType { STAT, SYST };

    std::pair<double, double> get_correlation(ParamId p1, ParamId p2);
    std::pair<double, double> get_correlation(Observables o1, Observables o2);
    double get_combined_correlation(ParamId p1, ParamId p2);
    double get_combined_correlation(Observables o1, Observables o2);
    
    void set_correlation_matrix(const SparseMatrix<ParamId>&& correlation_matrix, CorrelationType type);
    void set_correlation_matrix(const SparseMatrix<Observables>&& correlation_matrix, CorrelationType type);
    void merge_correlation_matrix(const SparseMatrix<ParamId>&& correlation_matrix, CorrelationType type);
    void merge_correlation_matrix(const SparseMatrix<Observables>&& correlation_matrix, CorrelationType type);

private:
    static SparseMatrix<ParamId> parameter_stat_correlations;
    static SparseMatrix<ParamId> parameter_syst_correlations;
    static SparseMatrix<Observables> observable_stat_correlations;
    static SparseMatrix<Observables> observable_syst_correlations;
};

#endif // __CORRELATIONREPO_H__
