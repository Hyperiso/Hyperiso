#ifndef KCHARGEDCURRENTSWILSONGROUP_H
#define KCHARGEDCURRENTSWILSONGROUP_H

#include "KChargedCurrentWilson.h"
#include "WilsonGroup.h"
#include "MartyWilson.h"


class KulnuCoefficientGroup : public CoefficientGroup {
public:
    KulnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<KulnuCoefficientGroup>(adapters, true); }

    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
        const std::unordered_map<std::string, std::shared_ptr<Block>>& src
    );
};

#endif