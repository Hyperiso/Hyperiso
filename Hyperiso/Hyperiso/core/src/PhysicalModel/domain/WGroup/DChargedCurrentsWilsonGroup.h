#ifndef DCHARGEDCURRENTSWILSONGROUP_H
#define DCHARGEDCURRENTSWILSONGROUP_H

#include "DChargedCurrentWilson.h"
#include "WilsonGroup.h"
#include "MartyWilson.h"


class DslnuCoefficientGroup : public CoefficientGroup {
public:
    DslnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<DslnuCoefficientGroup>(adapters, true); }

    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
        const BlockSrc& src
    );
};

class DdlnuCoefficientGroup : public CoefficientGroup {
public:
    DdlnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<DdlnuCoefficientGroup>(adapters, true); }

    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
        const BlockSrc& src
    );
};

#endif