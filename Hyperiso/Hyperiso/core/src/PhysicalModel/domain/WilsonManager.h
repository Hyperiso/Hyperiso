#if !defined(HYPERISO_WILSONMANAGER_H)
#define HYPERISO_WILSONMANAGER_H
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <set>
#include "WilsonSuper.h"
#include "WilsonGroupSuper.h"
#include "BWilsonGroupSuper.h"
#include "MemoryManager.h"
#include "QCDHelper.h"
#include "Utils.h"
#include "HyperisoMaster.h"
#include "ScaleSetter.h"

class CoefficientManager {
private:
    std::map<std::string, std::shared_ptr<CoefficientGroup>> coefficientGroups;

    bool has_bsm;
    std::string bsm_suffix;
    ParameterProxy wilson_p {ParameterType::WILSON};

    void throw_no_group_error(const std::string& groupName) const;

public:
    CoefficientManager() = default;

    void initialize(const std::string& lhaFile, Model model = Model::SM, 
                    bool use_marty = false, bool is_spectrum = false, bool has_wilsons = false, bool has_obs = false);

    static std::shared_ptr<CoefficientManager> Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order);

    void set_matching_scale(double mu_W);
    void set_hadronic_scale(double mu_h);

    void registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group);
    void init_group_matching(const std::string& groupName, const std::string& order);
    void init_group_hadronic(const std::string& groupName, const std::string& order);
    void switchbasis(const std::string& groupName);
    void update(std::string group, double mu_W, double mu_h);
    void post_init();
    void complete_wilson_block_with_op(WGroup group_id, ContributionType src_id, ContributionType dest_id, BWilsonBasis basis);
    void complete_wilson_block_from_copy(WGroup group_id, ContributionType src, ContributionType dest, BWilsonBasis basis=BWilsonBasis::STANDARD);
    
    void fill_sources_for_group(const std::string & groupName, const std::string& order, std::unordered_map<ParameterType, std::vector<std::string>>& src, int id=0);

    std::string getModel();
    std::shared_ptr<CoefficientGroup> getCoefficientGroup(const std::string& groupName) const;
    std::map<std::string, std::shared_ptr<CoefficientGroup>> getGroups();

    complex_t getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false);
    complex_t getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false);
    complex_t getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false);
    complex_t getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false);

    void printGroupCoefficients(const std::string& groupName) const;
};


#endif