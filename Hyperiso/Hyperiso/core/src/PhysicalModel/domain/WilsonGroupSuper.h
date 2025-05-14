#ifndef WILSON_GROUP_SUPER_H
#define WILSON_GROUP_SUPER_H

#include <unordered_map>
#include <optional>
#include "Include.h"
#include "Utils.h"
#include "WilsonSuper.h"
#include "BWilsonSuper.h"
#include "ChargedCurrentWilsonSuper.h"
#include "MartyWilsonSuper.h"
#include "ParameterProxy.h"
#include "UseMarty.h"
#include "BlockProxy.h"

using BRP = BWilsonRunningParameters;

class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:
    // Constructors
    CoefficientGroup() = default;
    CoefficientGroup(const CoefficientGroup&);
    CoefficientGroup(CoefficientGroup&&) = default;
    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs);

    void init(QCDOrder order);
    virtual void init_running_blocks(QCDOrder order) = 0;
    void init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, BWilsonBasis basis, bool inter, std::vector<ContributionType> type);

    virtual std::shared_ptr<CoefficientGroup> get_sm_group() {LOG_ERROR("LogicError", "cannot get_sm_group for non real group"); return nullptr;}
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

#endif