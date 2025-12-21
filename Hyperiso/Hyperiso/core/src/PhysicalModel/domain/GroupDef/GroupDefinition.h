#pragma once
#include <unordered_map>
#include <map>
#include <vector>
#include <functional>
#include "WilsonGroup.h"

enum class Backend { Builtin, Marty };

struct BuildContext {
    WilsonGroupAdapterConfig adapters;
    Model model   = Model::SM;
    Backend backend = Backend::Builtin;
    ContributionType contrib = ContributionType::SM;
    WGroupId group_id{};
    std::string group_name = "";
};

inline constexpr const char* MATCHING_BLOCK_PLACEHOLDER = "$MATCHING_BLOCK$";

struct GroupDefinition {
    WGroupId id;
    std::vector<WCoef> members;
    std::unordered_map<WilsonBasis, std::map<QCDOrder, CoefficientGroupSources>> sources;

    using SetupHook = std::function<void(const BuildContext&, CoefficientGroup&)>;
    std::unordered_map<Model, std::vector<SetupHook>> setup;
};

namespace GroupDefinitions {
    const GroupDefinition& B();
    const GroupDefinition& BPrime();
    const GroupDefinition& BScalar();
    const GroupDefinition& CC_bc();
    const GroupDefinition& CC_bu();
    const GroupDefinition& CC_cs();
    const GroupDefinition& CC_cd();
    const GroupDefinition& CC_su();
    const GroupDefinition& CC_du();
    const GroupDefinition& MesonMixing();
    const GroupDefinition& K();
    const GroupDefinition& get(WGroupId g);

    void register_custom(const GroupDefinition& def);
    bool has_custom(WGroupId id);
    const GroupDefinition& get_custom(WGroupId id);

    const GroupDefinition& get(WGroupId g);
}
