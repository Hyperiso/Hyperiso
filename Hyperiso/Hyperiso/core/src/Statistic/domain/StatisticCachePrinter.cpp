#include "StatisticCachePrinter.h"

#include <iostream>

void StatisticCachePrinter::print(const StatCache& cache, std::ostream& os) {
    for (const auto& elem : cache.eta_specs_real) {
        os << " eta_specs_real : " << elem.first << " = " << elem.second;
    }
    os << '\n';

    for (const auto& elem : cache.SigmaEta) {
        for (const auto& elem2 : elem.second) {
            os << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
        }
        os << '\n';
    }

    for (const auto& elem : cache.exp_obs) {
        os << " exp_obs : " << elem.first.str() << " = " << elem.second;
    }
    os << '\n';

    for (const auto& elem : cache.SigmaObs) {
        for (const auto& elem2 : elem.second) {
            os << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
        }
        os << '\n';
    }

    for (const auto& elem : cache.p_specs) {
        os << " p_specs : " << elem.first << " = " << elem.second;
    }
    os << '\n';

    os << "eta size : " << cache.eta_specs_real.size() << '\n';
    if (!cache.SigmaEta.empty()) {
        os << "etasigma size : " << cache.eta_specs_real.size()
           << " | " << cache.SigmaEta.begin()->second.size() << '\n';
    }
    os << "p_specs size : " << cache.p_specs.size() << '\n';
    os << "exp_obs : " << cache.exp_obs.size() << '\n';
    if (!cache.SigmaObs.empty()) {
        os << "SigmaObs size : " << cache.SigmaObs.size()
           << " | " << cache.SigmaObs.begin()->second.size() << '\n';
    }
}
