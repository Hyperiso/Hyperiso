#include "CorrelationRepo.h"

std::pair<double, double> CorrelationRepository::get_correlation(ParamId id1, ParamId id2) {
    auto key = std::make_pair(id1, id2);
    double stat = parameter_stat_correlations.contains(key) ? parameter_stat_correlations.at(key) : 0;
    double syst = parameter_syst_correlations.contains(key) ? parameter_syst_correlations.at(key) : 0;
    return {stat, syst};
}

std::pair<double, double> CorrelationRepository::get_correlation(Observables id1, Observables id2) {
    auto key = std::make_pair(id1, id2);
    double stat = observable_stat_correlations.contains(key) ? observable_stat_correlations.at(key) : 0;
    double syst = observable_syst_correlations.contains(key) ? observable_syst_correlations.at(key) : 0;
    return {stat, syst};
}

double CorrelationRepository::get_combined_correlation(ParamId p1, ParamId p2) {
    auto key = std::make_pair(p1, p2);
    double stat = parameter_stat_correlations.contains(key) ? parameter_stat_correlations.at(key) : 0;
    double syst = parameter_syst_correlations.contains(key) ? parameter_syst_correlations.at(key) : 0;
    return std::hypot(stat, syst);
}

double CorrelationRepository::get_combined_correlation(Observables o1, Observables o2) {
    auto key = std::make_pair(o1, o2);
    double stat = observable_stat_correlations.contains(key) ? observable_stat_correlations.at(key) : 0;
    double syst = observable_syst_correlations.contains(key) ? observable_syst_correlations.at(key) : 0;
    return std::hypot(stat, syst);
}

void CorrelationRepository::set_correlation_matrix(const SparseMatrix<ParamId> &&correlation_matrix, CorrelationType type) {
    if (type == CorrelationType::STAT) {
        parameter_stat_correlations = correlation_matrix;
    } else {
        parameter_syst_correlations = correlation_matrix;
    }   
}

void CorrelationRepository::set_correlation_matrix(const SparseMatrix<Observables> &&correlation_matrix, CorrelationType type) {
    if (type == CorrelationType::STAT) {
        observable_stat_correlations = correlation_matrix;
    } else {
        observable_syst_correlations = correlation_matrix;
    } 
}

void CorrelationRepository::merge_correlation_matrix(const SparseMatrix<ParamId> &&correlation_matrix, CorrelationType type) {
    for (const auto &corr: correlation_matrix) {
        if (type == CorrelationType::STAT) {
            parameter_stat_correlations.insert_or_assign(corr.first, corr.second);
        } else {
            parameter_syst_correlations.insert_or_assign(corr.first, corr.second);
        } 
    }
}

void CorrelationRepository::merge_correlation_matrix(const SparseMatrix<Observables> &&correlation_matrix, CorrelationType type) {
    for (const auto &corr: correlation_matrix) {
        if (type == CorrelationType::STAT) {
            observable_stat_correlations.insert_or_assign(corr.first, corr.second);
        } else {
            observable_syst_correlations.insert_or_assign(corr.first, corr.second);
        } 
    }
}
