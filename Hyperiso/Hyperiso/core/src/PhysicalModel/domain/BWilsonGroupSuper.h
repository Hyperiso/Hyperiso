#ifndef BWILSON_GROUP_SUPER_H
#define BWILSON_GROUP_SUPER_H

#include "BWilsonSuper.h"
#include "WilsonGroupSuper.h"
#include "MartyModelNameAPI.h"
#include "MartyModelPathAPI.h"


class BCoefficientGroup : public CoefficientGroup {
public:
    BCoefficientGroup();

    void set_gen(int new_gen) {}
    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BCoefficientGroup>(); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_2_LO_calculation   (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_1_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_2_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_1_NNLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);

    void init_running_parameter_blocks();
};


class BPrimeCoefficientGroup : public CoefficientGroup {
public:
    BPrimeCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BPrimeCoefficientGroup>(); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
};


class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BScalarCoefficientGroup>(); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_1_NLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
};


std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs);

#endif