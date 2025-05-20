#ifndef __CHARGEDCURRENTSWILSONGROUP_H__
#define __CHARGEDCURRENTSWILSONGROUP_H__

#include "ChargedCurrentWilson.h"
#include "WilsonGroupSuper.h"

class BlnuCoefficientGroup : public CoefficientGroup {
public:
    BlnuCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order);

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BlnuCoefficientGroup>(); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
};


class BclnuCoefficientGroup : public CoefficientGroup {
public:
    BclnuCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order);

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BclnuCoefficientGroup>(); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
};

#endif // __CHARGEDCURRENTSWILSONGROUP_H__
