#ifndef BWILSON_GROUP_SUPER_H
#define BWILSON_GROUP_SUPER_H

#include "BWilson.h"
#include "WilsonGroup.h"
#include "ICoreAPI.h"
#include "Block.h"
#include "BWilsonRunningParameters.h"
#include "MartyWilson.h"

using BRP = BWilsonRunningParameters;


class BCoefficientGroup : public CoefficientGroup {
public:
    BCoefficientGroup(WilsonGroupAdapterConfig adapters);

    void set_gen(int) {}
    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BCoefficientGroup>(adapters); }
    static std::unordered_map<WCoefId, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_2_LO_calculation   (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_1_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_2_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    static std::unordered_map<WCoefId, scalar_t> base_1_NNLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, const BlockSrc& src);
    
};




std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs);

#endif