#ifndef __CHARGEDCURRENTSWILSONGROUP_H__
#define __CHARGEDCURRENTSWILSONGROUP_H__

#include "ChargedCurrentWilsonSuper.h"
#include "WilsonGroupSuper.h"
#include "MartyModelPathAPI.h"
#include "MartyModelNameAPI.h"

class BlnuCoefficientGroup : public CoefficientGroup {
public:
    BlnuCoefficientGroup(bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BlnuCoefficientGroup>(true); }
    
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
        const std::unordered_map<std::string, std::shared_ptr<Block>>& src
    );

    void init_sources();
    void add_wilson_coefficients(bool force_sm=false);
};


class BclnuCoefficientGroup : public CoefficientGroup {
public:
    BclnuCoefficientGroup(bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BclnuCoefficientGroup>(true); }

private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
        const std::unordered_map<std::string, std::shared_ptr<Block>>& src
    );

    void init_sources();
    void add_wilson_coefficients(bool force_sm=false);
};

#endif // __CHARGEDCURRENTSWILSONGROUP_H__
