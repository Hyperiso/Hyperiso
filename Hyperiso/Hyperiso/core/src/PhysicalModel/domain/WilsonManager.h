#ifndef HYPERISO_WILSONMANAGER_H
#define HYPERISO_WILSONMANAGER_H
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <set>
#include "Wilson.h"
#include "WilsonGroup.h"
#include "BWilsonGroup.h"
#include "MemoryManager.h"
#include "QCDHelper.h"
#include "Utils.h"
#include "HyperisoMaster.h"
#include "IParamSetter.h"
#include "IWilsonParameters.h"

struct PortsConfig {

    PortsConfig(std::shared_ptr<IBlockComposer> iblock_c, std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy, std::shared_ptr<ICoreAPI<bool>> use_marty, std::shared_ptr<ICoreAPI<bool>> has_wilson, std::shared_ptr<ICoreAPI<Model>> model_api, std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api) :
        iblock_c(iblock_c), wilson_proxy(wilson_proxy), 
        use_marty(use_marty), has_wilson(has_wilson),
        model_api(model_api), scale_setter_api(scale_setter_api) {}

    std::shared_ptr<IBlockComposer> iblock_c;
    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy;
    std::shared_ptr<ICoreAPI<bool>> use_marty;
    std::shared_ptr<ICoreAPI<bool>> has_wilson;
    std::shared_ptr<ICoreAPI<Model>> model_api;
    std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api;
    std::function<std::shared_ptr<CoefficientGroup>(WGroupId, Model, bool, ContributionType, std::string)> build_group;
};

class CoefficientManager {
private:
    std::map<std::string, std::shared_ptr<CoefficientGroup>> coefficientGroups;
    PortsConfig ports_config;

    void throw_no_group_error(const std::string& groupName) const;

public:
    CoefficientManager(PortsConfig ports_config) : ports_config(ports_config) {}
    CoefficientManager(const CoefficientManager&) = delete;
    CoefficientManager operator=(const CoefficientManager&) = delete;

    static std::shared_ptr<CoefficientManager> Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order, PortsConfig portconfig, std::map<Model, std::shared_ptr<IWilsonParameterHelper>> wilson_param_helpers = {});

    void set_matching_scale(double mu_W);
    void set_hadronic_scale(double mu_h);

    void ensure_sm_intermediate_and_copy_to_final(
        const std::string& groupName,
        QCDOrder order,
        PortsConfig& ports_config
    );
    void compose_missing_from_calculation(
        const std::string& groupName,
        QCDOrder order,
        PortsConfig& ports_config
    );

    void compose_from_fwcoef(
        const std::string& groupName,
        QCDOrder order,
        PortsConfig& ports_config
    );
    void ensure_matching_triplet_zeroed(
        const std::string& groupName,
        QCDOrder o
    );
    void ensure_sm_model_triplet_in_matching(
        const std::string& groupName,
        QCDOrder max_order
    );
    void ensure_final_triplet_defaults_zero(
        const std::string& groupName,
        QCDOrder max_order
    );
    void registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group);
    void init_group_matching(const std::string& groupName, const std::string& order);
    void init_group_hadronic(const std::string& groupName, const std::string& order, WilsonBasis id);
    void init_group_hadronic_all_bases(const std::string& groupName, const std::string& order);
    void init_specific_order_group_matching(const std::string& groupName, const std::string& order, bool only_total);
    void update(std::string group, double mu_W, double mu_h);
    
    void fill_sources_for_group(const std::string & groupName, const std::string& order, std::unordered_map<ParameterType, std::vector<std::string>>& src, WilsonBasis id);
    void fill_matching_groups(const std::string& groupName, const std::string& order);

    std::string getModel();
    std::shared_ptr<CoefficientGroup> getCoefficientGroup(const std::string& groupName) const;
    std::map<std::string, std::shared_ptr<CoefficientGroup>> getGroups();
    std::unordered_set<WilsonBasis> getGroupBases(WGroupId group);

    complex_t getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type);
    complex_t getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type);
    complex_t getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type, WilsonBasis basis = WilsonBasis::B_STANDARD);
    complex_t getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type, WilsonBasis basis = WilsonBasis::B_STANDARD);

    void printGroupCoefficients(const std::string& groupName) const;

    ~CoefficientManager();
};


#endif