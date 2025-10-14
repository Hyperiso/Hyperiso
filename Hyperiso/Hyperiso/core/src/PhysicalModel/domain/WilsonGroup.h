#ifndef WILSON_GROUP_SUPER_H
#define WILSON_GROUP_SUPER_H

#include <unordered_map>
#include <optional>
#include "Include.h"
#include "Utils.h"
#include "Wilson.h"
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "IMartyWilsonProxy.h"
#include "config.hpp"
#include "InterpretedParam.h"

// #include "BWilson.h"
// #include "ChargedCurrentWilson.h"
// #include "MartyWilson.h"
// #include "ParameterProxy.h"
// #include "UseMarty.h"
// #include "BlockProxy.h"

struct CoefficientGroupSources {
    std::unordered_map<ParameterType, std::vector<std::string>> sources {};
    std::function<std::unordered_map<WCoef, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&, const std::unordered_map<std::string, std::shared_ptr<Block>>&)> func =
        [](const auto&, const auto&) { return std::unordered_map<WCoef, scalar_t>(); };
};

struct WilsonGroupAdapterConfig {

    WilsonGroupAdapterConfig(std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy, std::shared_ptr<IBlockComposer> ibc,
    std::shared_ptr<ICoreAPI<bool>> use_marty, std::shared_ptr<ICoreAPI<std::string>> marty_model_name, std::shared_ptr<ICoreAPI<fs::path>> marty_model_path,
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr) 
    : wilson_proxy(wilson_proxy), iblock_c(ibc), use_marty(use_marty), marty_model_name(marty_model_name), marty_model_path(marty_model_path), marty_proxy(marty_proxy) {}

    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy;
    std::shared_ptr<IBlockComposer> iblock_c;
    std::shared_ptr<ICoreAPI<bool>> use_marty;
    std::shared_ptr<ICoreAPI<std::string>> marty_model_name;
    std::shared_ptr<ICoreAPI<fs::path>> marty_model_path;
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;

    fs::path sm_path = fs::path(std::string(project_assets_root.data())+"input_files/marty_model/sm.h");
};

class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:
    // Constructors
    CoefficientGroup(WilsonGroupAdapterConfig adapters) : adapters(adapters) {};
    CoefficientGroup(const CoefficientGroup&);
    CoefficientGroup(CoefficientGroup&&) = default;
    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs, WilsonGroupAdapterConfig adapters);

    void init(QCDOrder order);
    void init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, WilsonBasis basis, bool inter, std::vector<ContributionType> type);

    virtual std::shared_ptr<CoefficientGroup> get_sm_group() {LOG_ERROR("LogicError", "cannot get_sm_group for virtual group"); return nullptr;}

    // Getters
    complex_t get_matching_coefficient(std::string coeff, std::string order, ContributionType cont_type) const;
    complex_t get_running_coefficient(std::string coeff, std::string order, ContributionType cont_type, WilsonBasis basis = WilsonBasis::B_STANDARD) const;
    QCDOrder get_order();
    std::string get_matching_storage_block() const { return GroupMapper::str(this->id, ScaleType::MATCHING); }
    ContributionType get_type() {return this->wilson_type;}
    std::unordered_map<ParameterType, std::vector<std::string>> get_sources(QCDOrder ord, WilsonBasis id) {return this->sources[id][ord].sources;}
    std::unordered_set<WilsonBasis> get_bases() const {
        return get_keys(this->sources);
    }
    std::function<std::unordered_map<WCoef, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&, const std::unordered_map<std::string, std::shared_ptr<Block>>&)> get_func(QCDOrder ord, WilsonBasis id) {return this->sources[id][ord].func;}

    // Interface methods
    virtual std::shared_ptr<CoefficientGroup> clone() const = 0;

    virtual ~CoefficientGroup() = default;

protected:
    void claim_coefficients();

    // static complex_t ensure_coef(WCoef coef, QCDOrder order, ContributionType type, std::string matching_block);
    ContributionType wilson_type {ContributionType::SM};
    QCDOrder current_order = QCDOrder::LO;
    WGroup id;
    std::map<WilsonBasis, std::map<QCDOrder, CoefficientGroupSources>> sources;

    WilsonGroupAdapterConfig adapters;
};

std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);
std::ostream& operator<<(std::ostream& os, const std::shared_ptr<CoefficientGroup>& coeffs);
#endif