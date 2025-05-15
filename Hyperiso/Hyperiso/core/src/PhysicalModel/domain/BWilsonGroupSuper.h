#ifndef BWILSON_GROUP_SUPER_H
#define BWILSON_GROUP_SUPER_H

#include "BWilsonSuper.h"
#include "WilsonGroupSuper.h"

class BCoefficientGroup : public CoefficientGroup {
public:
    BCoefficientGroup();

    void init_running_blocks(QCDOrder order) override;
    void set_gen(int new_gen) {}
    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BCoefficientGroup>(); }
private:
    static void base_1_LO_calculation   (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
    static void base_2_LO_calculation   (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
    static void base_1_NLO_calculation  (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
    static void base_2_NLO_calculation  (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
    static void base_1_NNLO_calculation (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);

    void init_running_parameter_blocks();
};


class BPrimeCoefficientGroup : public CoefficientGroup {
public:
    BPrimeCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order) override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BPrimeCoefficientGroup>(); }
private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order) override;

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BScalarCoefficientGroup>(); }
private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
    static void base_1_NLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


class BlnuCoefficientGroup : public CoefficientGroup {
public:
    BlnuCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order);

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BlnuCoefficientGroup>(); }
private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


class BclnuCoefficientGroup : public CoefficientGroup {
public:
    BclnuCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order);

    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BclnuCoefficientGroup>(); }
private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs);

#endif