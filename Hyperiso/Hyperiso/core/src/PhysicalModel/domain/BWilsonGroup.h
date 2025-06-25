#ifndef BWILSON_GROUP_SUPER_H
#define BWILSON_GROUP_SUPER_H

#include "BWilson.h"
#include "WilsonGroup.h"
#include "MartyModelNameAPI.h"
#include "MartyModelPathAPI.h"


class BCoefficientGroup : public CoefficientGroup {
public:
    BCoefficientGroup(bool force_sm=false);

    void set_gen(int new_gen) {}
    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BCoefficientGroup>(true); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_2_LO_calculation   (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_1_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_2_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_1_NNLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);

    void init_running_parameter_blocks();

    void init_sources();
    void add_wilson_coefficients(bool force_sm=false);
};


class BPrimeCoefficientGroup : public CoefficientGroup {
public:
    BPrimeCoefficientGroup(bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BPrimeCoefficientGroup>(true); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);

    void init_sources();
    void add_wilson_coefficients(bool force_sm=false);
};


class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup(bool force_sm=false);
    std::shared_ptr<CoefficientGroup> clone() const override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BScalarCoefficientGroup>(true); }
private:
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    static std::unordered_map<WCoef, scalar_t> base_1_NLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);

    void init_sources();
    void add_wilson_coefficients(bool force_sm=false);
};


std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs);

#endif