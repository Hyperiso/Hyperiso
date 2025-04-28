#include "CorrelationProvider.h"

double CorrelationProvider::operator()(const ParamId &pid_1, const ParamId &pid_2, CorrelationType type) {
    return get_correlation(pid_1, pid_2, type);
}

double CorrelationProvider::operator()(const Observables &oid_1, const Observables &oid_2, CorrelationType type) {
    return get_correlation(oid_1, oid_2, type);
}

template <typename T>
inline double CorrelationProvider::get_correlation(const T &id_1, const T &id_2, CorrelationType type) const {
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