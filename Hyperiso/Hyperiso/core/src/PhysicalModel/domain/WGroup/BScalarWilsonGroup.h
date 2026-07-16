#ifndef BSCALAR_WILSON_GROUP_H
#define BSCALAR_WILSON_GROUP_H

#include "WilsonGroup.h"

class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup(WilsonGroupAdapterConfig adapters);
    std::shared_ptr<CoefficientGroup> clone() const override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BScalarCoefficientGroup>(adapters); }
    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_1_NLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    
};

#endif