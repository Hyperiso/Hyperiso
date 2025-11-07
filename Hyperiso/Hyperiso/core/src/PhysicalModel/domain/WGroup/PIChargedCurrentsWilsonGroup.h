#ifndef PICHARGEDCURRENTSWILSONGROUP_H
#define PICHARGEDCURRENTSWILSONGROUP_H

#include "PIChargedCurrentWilson.h"
#include "WilsonGroup.h"
#include "MartyWilson.h"


class PIulnuCoefficientGroup : public CoefficientGroup {
public:
    PIulnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<PIulnuCoefficientGroup>(adapters, true); }

    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& src
    );
};

#endif