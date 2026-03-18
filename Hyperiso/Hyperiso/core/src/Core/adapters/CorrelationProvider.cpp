#include "CorrelationProvider.h"

double CorrelationProvider::operator()(const ParamId& pid_1,
                                       const ParamId& pid_2,
                                       CorrelationType type) {
    return get_correlation(pid_1, pid_2, type);
}

double CorrelationProvider::operator()(const ExperimentObs& oid_1,
                                       const ExperimentObs& oid_2,
                                       CorrelationType type) {
    return get_correlation(oid_1, oid_2, type);
}

double CorrelationProvider::operator()(const std::string& experiment,
                                       const Observables& oid_1,
                                       const Observables& oid_2,
                                       CorrelationType type) {
    ObservableId obs_id_1 = ObservableMapper::to_id(oid_1);
    ObservableId obs_id_2 = ObservableMapper::to_id(oid_2);
    return (*this)(experiment, obs_id_1, obs_id_2, type);
}

double CorrelationProvider::operator()(const std::string& experiment,
                                       const ObservableId& oid_1,
                                       const ObservableId& oid_2,
                                       CorrelationType type) {
    BinnedObservableId oid_1_b{oid_1, {0., 0.}};
    BinnedObservableId oid_2_b{oid_2, {0., 0.}};
    return (*this)(experiment, oid_1_b, oid_2_b, type);
}

double CorrelationProvider::operator()(const std::string& experiment,
                                       const BinnedObservableId& oid_1,
                                       const BinnedObservableId& oid_2,
                                       CorrelationType type) {
    return get_correlation(experiment, oid_1, oid_2, type);
}

double CorrelationProvider::operator()(const std::string& exp_1,
                                       const BinnedObservableId& oid_1,
                                       const std::string& exp_2,
                                       const BinnedObservableId& oid_2,
                                       CorrelationType type) {
    return get_correlation(exp_1, oid_1, exp_2, oid_2, type);
}

bool CorrelationProvider::exists(const ParamId& pid_1,
                                 const ParamId& pid_2,
                                 CorrelationType type) {
    return get_correlation(pid_1, pid_2, type) != 0.0;
}

bool CorrelationProvider::exists(const ExperimentObs& oid_1,
                                 const ExperimentObs& oid_2,
                                 CorrelationType type) {
    return get_correlation(oid_1, oid_2, type) != 0.0;
}

bool CorrelationProvider::exists(const std::string& experiment,
                                 const Observables& oid_1,
                                 const Observables& oid_2,
                                 CorrelationType type) {
    return (*this)(experiment, oid_1, oid_2, type) != 0.0;
}

bool CorrelationProvider::exists(const std::string& experiment,
                                 const ObservableId& oid_1,
                                 const ObservableId& oid_2,
                                 CorrelationType type) {
    return (*this)(experiment, oid_1, oid_2, type) != 0.0;
}

bool CorrelationProvider::exists(const std::string& experiment,
                                 const BinnedObservableId& oid_1,
                                 const BinnedObservableId& oid_2,
                                 CorrelationType type) {
    return (*this)(experiment, oid_1, oid_2, type) != 0.0;
}

bool CorrelationProvider::exists(const std::string& exp_1,
                                 const BinnedObservableId& oid_1,
                                 const std::string& exp_2,
                                 const BinnedObservableId& oid_2,
                                 CorrelationType type) {
    return (*this)(exp_1, oid_1, exp_2, oid_2, type) != 0.0;
}

template <typename T>
inline double CorrelationProvider::get_correlation(const T& id_1,
                                                   const T& id_2,
                                                   CorrelationType type) const {
    if (id_1 == id_2) {
        return 1.0;
    }

    CorrelationRepository correlation_repo = MemoryManager::GetInstance()->get_correlation_repository();

    switch (type) {
    case CorrelationType::STAT:
        return correlation_repo.get_correlation(id_1, id_2).first;
    case CorrelationType::SYST:
        return correlation_repo.get_correlation(id_1, id_2).second;
    case CorrelationType::COMBINED:
        return correlation_repo.get_combined_correlation(id_1, id_2);
    default:
        LOG_ERROR("CorrelationProvider", "Unknown correlation type requested.");
        return 0.0;
    }
}

double CorrelationProvider::get_correlation(const std::string& experiment,
                                            const BinnedObservableId& id_1,
                                            const BinnedObservableId& id_2,
                                            CorrelationType type) const {
    if (id_1 == id_2) {
        return 1.0;
    }

    CorrelationRepository correlation_repo = MemoryManager::GetInstance()->get_correlation_repository();

    switch (type) {
    case CorrelationType::STAT:
        return correlation_repo.get_correlation(experiment, id_1, id_2).first;
    case CorrelationType::SYST:
        return correlation_repo.get_correlation(experiment, id_1, id_2).second;
    case CorrelationType::COMBINED: {
        auto corr = correlation_repo.get_correlation(experiment, id_1, id_2);
        return std::hypot(corr.first, corr.second);
    }
    default:
        LOG_ERROR("CorrelationProvider", "Unknown correlation type requested.");
        return 0.0;
    }
}

double CorrelationProvider::get_correlation(const std::string& exp_1,
                                            const BinnedObservableId& id_1,
                                            const std::string& exp_2,
                                            const BinnedObservableId& id_2,
                                            CorrelationType type) const {
    if (exp_1 == exp_2 && id_1 == id_2) {
        return 1.0;
    }

    CorrelationRepository correlation_repo = MemoryManager::GetInstance()->get_correlation_repository();

    switch (type) {
    case CorrelationType::STAT:
        return correlation_repo.get_correlation(exp_1, id_1, exp_2, id_2).first;
    case CorrelationType::SYST:
        return correlation_repo.get_correlation(exp_1, id_1, exp_2, id_2).second;
    case CorrelationType::COMBINED: {
        auto corr = correlation_repo.get_correlation(exp_1, id_1, exp_2, id_2);
        return std::hypot(corr.first, corr.second);
    }
    default:
        LOG_ERROR("CorrelationProvider", "Unknown correlation type requested.");
        return 0.0;
    }
}