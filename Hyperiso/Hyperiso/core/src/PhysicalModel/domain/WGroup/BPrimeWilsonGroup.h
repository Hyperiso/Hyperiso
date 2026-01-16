#ifndef BPRIME_WILSON_GROUP_H
#define BPRIME_WILSON_GROUP_H

#include "WilsonGroup.h"
#include "BWilsonRunningParameters.h"

using BRP = BWilsonRunningParameters;

class BPrimeCoefficientGroup : public CoefficientGroup {
public:
    BPrimeCoefficientGroup(WilsonGroupAdapterConfig adapters);
    std::shared_ptr<CoefficientGroup> clone() const override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BPrimeCoefficientGroup>(adapters); }
    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);

};

#endif