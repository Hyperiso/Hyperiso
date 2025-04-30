#ifndef __CONFIGS_H__
#define __CONFIGS_H__

#include "Include.h"

struct WilsonBuildConfig : public AbstractConfig {
    std::unordered_set<WGroup> groups;
    double matching_scale;
    double hadronic_scale;
    QCDOrder order {QCDOrder::LO}; 
};

struct WilsonRequest : public AbstractConfig {
    WGroup group;
    WCoef coefficient;
    QCDOrder order {QCDOrder::LO};
    ContributionType contribution {ContributionType::TOTAL};
    ScaleType scale_type {ScaleType::HADRONIC};
    bool sum_qcd_orders {false};
};

struct AlphasConfig : public AbstractConfig {
    double scale;
    MassType m_b_type;
    MassType m_t_type;

    AlphasConfig(double scale, MassType m_b_type, MassType m_t_type) : scale(scale), m_b_type(m_b_type), m_t_type(m_t_type) {}
};

struct MassConfig : public AlphasConfig {
    int pdg_id;

    MassConfig(int pdg_id, double scale, MassType m_b_type, MassType m_t_type) : AlphasConfig(scale, m_b_type, m_t_type), pdg_id(pdg_id) {}
};

#endif // __CONFIGS_H__
