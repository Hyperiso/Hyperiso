#ifndef WILSON_GROUP_H
#define WILSON_GROUP_H

#include <unordered_map>
#include <optional>
#include "Include.h"
#include "Utils.h"
#include "Wilson.h"
#include "BWilson.h"
#include "ChargedCurrentWilson.h"
#include "MartyWilson.h"
#include "ParameterProxy.h"
#include "UseMarty.h"

using BRP = BWilsonRunningParameters;

class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:
    // Constructors
    CoefficientGroup() = default;
    CoefficientGroup(const CoefficientGroup&) = default;
    CoefficientGroup(CoefficientGroup&&) = default;
    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs);

    void init(QCDOrder order);
    virtual void init_running_blocks(QCDOrder order) = 0;
    void init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, BWilsonBasis basis, bool inter, std::vector<ContributionType> type);

    // Getters
    complex_t get_matching_coefficient(std::string coeff, std::string order) const;
    complex_t get_running_coefficient(std::string coeff, std::string order) const;
    QCDOrder get_order();
    bool is_double_basis() const;

    // Interface methods
    virtual std::shared_ptr<CoefficientGroup> clone() const = 0;
    void switch_basis();

    virtual ~CoefficientGroup() = default;

protected:
    void claim_coefficients();

    static complex_t ensure_coef(WCoef coef, QCDOrder order, ContributionType type, std::string matching_block);

    std::optional<BWilsonBasis> basis;
    ContributionType wilson_type {ContributionType::SM};
    QCDOrder current_order = QCDOrder::LO; //TODO SAME : cannot be none, need to see logic
    WGroup id;
};


class BCoefficientGroup : public CoefficientGroup {
public:
    BCoefficientGroup();

    void init_running_blocks(QCDOrder order) override;
    void set_gen(int new_gen) {}
    std::shared_ptr<CoefficientGroup> clone() const override;

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

private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order) override;

private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
    static void base_1_NLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


class BlnuCoefficientGroup : public CoefficientGroup {
public:
    BlnuCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order);

private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


class BclnuCoefficientGroup : public CoefficientGroup {
public:
    BclnuCoefficientGroup();
    std::shared_ptr<CoefficientGroup> clone() const override;
    void init_running_blocks(QCDOrder order);

private:
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
};


std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs);

#endif