#ifndef __MESONMIXINGWILSONGROUP_H__
#define __MESONMIXINGWILSONGROUP_H__

#include "MesonMixingWilson.h"
#include "WilsonGroup.h"
// #include "MartyModelNameAPI.h"
// #include "MartyModelPathAPI.h"
#include "MesonMixingRunningParameters.h"
#include "MartyWilson.h"
#include "QCDHelper.h"


class MesonMixingCoefficientGroup : public CoefficientGroup {
public:
    MesonMixingCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm=false);

    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<MesonMixingCoefficientGroup>(adapters, true); }
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    
private:

    void init_running_parameter_blocks();

    // void init_sources();
    // void add_wilson_coefficients(bool force_sm=false);
};

#endif // __MESONMIXINGWILSONGROUP_H__
