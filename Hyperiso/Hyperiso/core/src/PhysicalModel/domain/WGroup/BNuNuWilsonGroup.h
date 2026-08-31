#ifndef BNUNU_WILSON_GROUP_H
#define BNUNU_WILSON_GROUP_H

#include "WilsonGroup.h"

class BNuNuCoefficientGroup : public CoefficientGroup {
public:
    BNuNuCoefficientGroup(WilsonGroupAdapterConfig adapters);

    void set_gen(int) {}
    std::shared_ptr<CoefficientGroup> clone() const override;
    std::shared_ptr<CoefficientGroup> get_sm_group() override {
        return std::make_shared<BNuNuCoefficientGroup>(adapters);
    }

    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_1_NLO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_1_NNLO_calculation(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& src);

private:
    void init_sources() {}
    void add_wilson_coefficients() {}
};

#endif
