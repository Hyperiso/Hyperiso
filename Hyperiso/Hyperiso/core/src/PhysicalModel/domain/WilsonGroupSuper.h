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

struct CoefficientGroupSources {
    std::unordered_map<ParameterType, std::vector<std::string>> sources {};
    std::function<std::unordered_map<WCoef, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&, const std::unordered_map<std::string, std::shared_ptr<Block>>&)> func =
        [](const auto&, const auto&) { return std::unordered_map<WCoef, scalar_t>(); };
};


class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:
    // Constructors
    CoefficientGroup() = default;
    CoefficientGroup(const CoefficientGroup&);
    CoefficientGroup(CoefficientGroup&&) = default;
    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs);

    void init(QCDOrder order);
    void init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, BWilsonBasis basis, bool inter, std::vector<ContributionType> type);

    virtual std::shared_ptr<CoefficientGroup> get_sm_group() {LOG_ERROR("LogicError", "cannot get_sm_group for virtual group"); return nullptr;}

    // Getters
    complex_t get_matching_coefficient(std::string coeff, std::string order, ContributionType cont_type) const;
    complex_t get_running_coefficient(std::string coeff, std::string order, ContributionType cont_type) const;
    QCDOrder get_order();
    std::string get_matching_storage_block() const { return GroupMapper::str(this->id, ScaleType::MATCHING); }
    ContributionType get_type() {return this->wilson_type;}
    std::unordered_map<ParameterType, std::vector<std::string>> get_sources(QCDOrder ord, int id=0) {return this->sources[id][ord].sources;}
    std::function<std::unordered_map<WCoef, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&, const std::unordered_map<std::string, std::shared_ptr<Block>>&)> get_func(QCDOrder ord, int id=0) {return this->sources[id][ord].func;}
    
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
    QCDOrder current_order = QCDOrder::LO;
    WGroup id;
    std::vector<std::map<QCDOrder, CoefficientGroupSources>> sources;
};

#endif