#ifndef BCHARGEDCURRENTSWILSONGROUP_H
#define BCHARGEDCURRENTSWILSONGROUP_H

#include "ChargedCurrentWilson.h"
#include "WilsonGroup.h"
#include "MartyWilson.h"


class BclnuCoefficientGroup : public CoefficientGroup {
public:
    BclnuCoefficientGroup(WilsonGroupAdapterConfig adapters);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BclnuCoefficientGroup>(adapters); }

    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& src
    );
};

class BulnuCoefficientGroup : public CoefficientGroup {
public:
    BulnuCoefficientGroup(WilsonGroupAdapterConfig adapters);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BulnuCoefficientGroup>(adapters); }

    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& src
    );
};

#endif
