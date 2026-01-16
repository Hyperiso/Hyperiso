#ifndef MESONMIXINGWILSONGROUP_H
#define MESONMIXINGWILSONGROUP_H

#include "MesonMixingWilson.h"
#include "WilsonGroup.h"
#include "MesonMixingRunningParameters.h"
#include "MartyWilson.h"
#include "QCDHelper.h"


class MesonMixingCoefficientGroup : public CoefficientGroup {
public:
    MesonMixingCoefficientGroup(WilsonGroupAdapterConfig adapters);

    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<MesonMixingCoefficientGroup>(adapters); }
    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    
private:
    void init_running_parameter_blocks();
};

#endif // MESONMIXINGWILSONGROUP_H
