// #ifndef WILSON_GROUP_H
// #define WILSON_GROUP_H

// #include <unordered_map>
// #include <optional>
// #include "Include.h"
// #include "Utils.h"
// #include "Wilson.h"
// #include "BWilson.h"
// #include "ChargedCurrentWilson.h"
// #include "MartyWilson.h"
// #include "ParameterProxy.h"
// #include "UseMarty.h"
// #include "BlockProxy.h"

// using BRP = BWilsonRunningParameters;

// class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
// public:
//     // Constructors
//     CoefficientGroup() = default;
//     CoefficientGroup(const CoefficientGroup&);
//     CoefficientGroup(CoefficientGroup&&) = default;
//     CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs);

//     void init(QCDOrder order);
//     virtual void init_running_blocks(QCDOrder order) = 0;
//     void init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, WilsonBasis basis, bool inter, std::vector<ContributionType> type);

//     virtual std::shared_ptr<CoefficientGroup> get_sm_group() {LOG_ERROR("LogicError", "cannot get_sm_group for non real group"); return nullptr;}
//     // Getters
//     complex_t get_matching_coefficient(std::string coeff, std::string order) const;
//     complex_t get_running_coefficient(std::string coeff, std::string order) const;
//     QCDOrder get_order();
//     bool is_double_basis() const;

//     // Interface methods
//     virtual std::shared_ptr<CoefficientGroup> clone() const = 0;
//     void switch_basis();

//     virtual ~CoefficientGroup() = default;

// protected:
//     void claim_coefficients();

//     static complex_t ensure_coef(WCoef coef, QCDOrder order, ContributionType type, std::string matching_block);

//     std::optional<WilsonBasis> basis;
//     ContributionType wilson_type {ContributionType::SM};
//     QCDOrder current_order = QCDOrder::LO; //TODO SAME : cannot be none, need to see logic
//     WGroup id;
// };


// class BCoefficientGroup : public CoefficientGroup {
// public:
//     BCoefficientGroup();

//     void init_running_blocks(QCDOrder order) override;
//     void set_gen(int new_gen) {}
//     std::shared_ptr<CoefficientGroup> clone() const override;
    
//     std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BCoefficientGroup>(); }
// private:
//     static void base_1_LO_calculation   (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
//     static void base_2_LO_calculation   (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
//     static void base_1_NLO_calculation  (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
//     static void base_2_NLO_calculation  (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
//     static void base_1_NNLO_calculation (const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);

//     void init_running_parameter_blocks();
// };


// class BPrimeCoefficientGroup : public CoefficientGroup {
// public:
//     BPrimeCoefficientGroup();
//     std::shared_ptr<CoefficientGroup> clone() const override;
//     void init_running_blocks(QCDOrder order) override;

//     std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BPrimeCoefficientGroup>(); }
// private:
//     static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
// };


// class BScalarCoefficientGroup : public CoefficientGroup {
// public:
//     BScalarCoefficientGroup();
//     std::shared_ptr<CoefficientGroup> clone() const override;
//     void init_running_blocks(QCDOrder order) override;

//     std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BScalarCoefficientGroup>(); }
// private:
//     static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
//     static void base_1_NLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
// };


// class BlnuCoefficientGroup : public CoefficientGroup {
// public:
//     BlnuCoefficientGroup();
//     std::shared_ptr<CoefficientGroup> clone() const override;
//     void init_running_blocks(QCDOrder order);

//     std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BlnuCoefficientGroup>(); }
// private:
//     static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
// };


// class BclnuCoefficientGroup : public CoefficientGroup {
// public:
//     BclnuCoefficientGroup();
//     std::shared_ptr<CoefficientGroup> clone() const override;
//     void init_running_blocks(QCDOrder order);

//     std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<BclnuCoefficientGroup>(); }
// private:
//     static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>, ContributionType);
// };


// std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
// std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs);

// #endif