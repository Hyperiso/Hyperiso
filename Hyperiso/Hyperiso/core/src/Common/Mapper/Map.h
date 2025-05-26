#ifndef MAP_H
#define MAP_H

#include "GeneralEnum.h"
#include "ObservableMapper.h"
#include "QCDOrderMapper.h"
#include "GroupMapper.h"
#include "WCoeffMapper.h"
#include "ParameterTypeMapper.h"
#include "ModelMapper.h"
#include "WilsonBasisMapper.h"
#include "ContributionTypeMapper.h"
#include "MassTypeMapper.h"
#include "ScaleTypeMapper.h"
#include "DecayMapper.h"

const std::map<Observables, std::string>& observable_mapping();
const std::map<QCDOrder, std::string>& order_mapping();
const std::map<WGroup, std::string>& group_mapping();
const std::map<WCoef, std::string>& wcoef_mapping();
const std::map<WCoef, std::pair<int, int>>& wcoef_flha_mapping();
const std::map<ParameterType, std::string>& parametertype_mapping();
const std::map<Model, std::string>& model_mapping();
const std::map<WilsonBasis, std::string>& wilsonbasis_mapping();
const std::map<ContributionType, std::string>& contributiontype_mapping();
const std::map<Decays, std::string>& decays_mapping();

const std::map<Decays, std::vector<Observables>>& decay_observable_mapping();

// inline const std::map<Observables, std::string>& ObservableMapper::mapping() {
//     static const std::map<Observables, std::string> m = observable_mapping();
//     return m;
// }

// inline const std::map<QCDOrder, std::string>& OrderMapper::mapping() {
//     static const std::map<QCDOrder, std::string> m = order_mapping();
//     return m;
// }

// inline const std::map<WGroup, std::string>& GroupMapper::mapping() {
//     static const std::map<WGroup, std::string> m = group_mapping();
//     return m;
// }

inline const std::map<WCoef, std::string>& WCoefMapper::_mapping() {
    static const std::map<WCoef, std::string> m = wcoef_mapping();
    return m;
}

inline const std::map<WCoef, std::pair<int, int>>& WCoefMapper::_flha_mapping() {
    static const std::map<WCoef, std::pair<int, int>> m = wcoef_flha_mapping();
    return m;
}

inline const std::map<ParameterType, std::string>& ParameterTypeMapper::mapping() {
    static const std::map<ParameterType, std::string> m = parametertype_mapping();
    return m;
}

inline const std::map<std::string, ParameterType>& ParameterTypeMapper::inverse_mapping() {
    static const std::map<std::string, ParameterType> inv = invert_map(mapping());
    return inv;
}

inline const std::map<Model, std::string>& ModelMapper::mapping() {
    static const std::map<Model, std::string> m = model_mapping();
    return m;
}

inline const std::map<std::string, Model>& ModelMapper::inverse_mapping() {
    static const std::map<std::string, Model> inv = invert_map(mapping());
    return inv;
}

inline const std::map<WilsonBasis, std::string>& WilsonBasisMapper::mapping() {
    static const std::map<WilsonBasis, std::string> m = {
        {WilsonBasis::B_STANDARD, "STANDARD"},
        {WilsonBasis::B_TRADITIONAL, "TRADITIONAL"},
    };
    return m;
}

inline const std::map<std::string, WilsonBasis>& WilsonBasisMapper::inverse_mapping() {
    static const std::map<std::string, WilsonBasis> inv = invert_map(mapping());
    return inv;
}

inline const std::map<ContributionType, std::string>& ContributionTypeMapper::mapping() {
    static const std::map<ContributionType, std::string> m = {
        {ContributionType::SM, "SM"},
        {ContributionType::BSM, "BSM"},
        {ContributionType::TOTAL, "TOTAL"},
    };
    return m;
}

inline const std::map<std::string, ContributionType>& ContributionTypeMapper::inverse_mapping() {
    static const std::map<std::string, ContributionType> inv = invert_map(mapping());
    return inv;
}

inline const std::map<MassType, std::string>& MassTypeMapper::mapping() {
    static const std::map<MassType, std::string> m = {
        {MassType::POLE, "POLE"},
        {MassType::MSBAR, "MSBAR"},
    };
    return m;
}

inline const std::map<std::string, MassType>& MassTypeMapper::inverse_mapping() {
    static const std::map<std::string, MassType> inv = invert_map(mapping());
    return inv;
}

inline const std::map<ScaleType, std::string>& ScaleTypeMapper::mapping() {
    static const std::map<ScaleType, std::string> m = {
        {ScaleType::MATCHING, "EW_SCALE"},
        {ScaleType::HADRONIC, "B_SCALE"},
    };
    return m;
}

inline const std::map<std::string, ScaleType>& ScaleTypeMapper::inverse_mapping() {
    static const std::map<std::string, ScaleType> inv = invert_map(mapping());
    return inv;
}

inline const std::map<Decays, std::string>& DecayMapper::mapping() {
    static const std::map<Decays, std::string> m = decays_mapping();
    return m;
}

inline const std::map<std::string, Decays>& DecayMapper::inverse_mapping() {
    static const std::map<std::string, Decays> inv = invert_map(mapping());
    return inv;
}

inline std::vector<Observables> DecayMapper::get_observables(Decays decay) {
    return decay_observable_mapping().at(decay);
}

inline Decays DecayMapper::get_decay(Observables obs) {
    for (const auto& [k, v] : decay_observable_mapping()) {
        if (std::find(v.begin(), v.end(), obs) != v.end()) {
            return k;
        }
    }
    LOG_ERROR("DecayMapper", "Observable belongs to none of the implemented decays.");
}

#endif