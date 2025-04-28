#include "CorrelationRepo.h"

std::pair<double, double> CorrelationRepository::get_correlation(ParamId id1, ParamId id2) const {
    auto key = std::make_pair(id1, id2);
    return parameter_correlations->at(key);
}

std::pair<double, double> CorrelationRepository::get_correlation(Observables id1, Observables id2) const {
    auto key = std::make_pair(id1, id2);
    return observable_correlations->at(key);
}

void CorrelationRepository::set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrices) {
    parameter_correlations = correlation_matrices;
}

void CorrelationRepository::set_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<Observables>> correlation_matrices) {
    observable_correlations = correlation_matrices;
}

void CorrelationRepository::merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<ParamId>> correlation_matrices) {
    for (const auto &corr: correlation_matrices->stat) {
        parameter_correlations->stat.insert_or_assign(corr.first, corr.second);
    }

    for (const auto &corr: correlation_matrices->syst) {
        parameter_correlations->syst.insert_or_assign(corr.first, corr.second);
    }
}

void CorrelationRepository::merge_correlation_matrix(std::shared_ptr<CorrelationMatrixPair<Observables>> correlation_matrices) {
    for (const auto &corr: correlation_matrices->stat) {
        observable_correlations->stat.insert_or_assign(corr.first, corr.second);
    }

    for (const auto &corr: correlation_matrices->syst) {
        observable_correlations->syst.insert_or_assign(corr.first, corr.second);
    }
}

void CorrelationRepository::print_content() const {
    LOG_INFO("------- Parameter correlations (off-diagonal) -------");
    LOG_INFO(parameter_correlations);
    LOG_INFO("------- Observable correlations (off-diagonal) -------");
    LOG_INFO(observable_correlations);
}

template<typename U>
std::ostream &operator<<(std::ostream &os, const CorrelationMatrixPair<U>& cmp) {
    for (const std::pair<std::pair<U, U>, double>& corr : cmp.stat) {
        auto corr_vals = cmp.at(corr.first);
        os << "\t(" << corr.first.first << ", " << corr.first.second << "): " << corr_vals.first << " (stat) + " << corr_vals.second << " (syst)\n";
    }
    return os;
}
