#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "WilsonInterface.h"
#include "mapper_hub.hpp"

namespace {

// Helper dev: resolve Wilson groups using runtime names instead of hard-coding WGroup enums.
// Important: GroupMapper uses the canonical names from Common/Mapper/Map.cpp, e.g.
// "BCoefficients", "BPrimeCoefficients", "BScalarCoefficients". You can add shorter
// aliases for your own custom groups with GroupMapper::register_custom(...).
std::unordered_set<WGroupId> resolve_groups(const std::vector<std::string>& names) {
    std::unordered_set<WGroupId> groups;

    for (const auto& name : names) {
        WGroupId gid = GroupMapper::id_of(name);
        std::cout << "Resolved Wilson group: " << name << " -> " << GroupMapper::str(gid) << "\n";
        groups.insert(gid);
    }

    return groups;
}

// Helper dev: query a Wilson coefficient by runtime ids.
// This uses the dynamic overloads added to WilsonInterface: no conversion back to WGroup/WCoef
// is required, so the same pattern also works for custom groups/coefficients.
scalar_t get_full_running_by_name(WilsonInterface& wilson,
                                  const std::string& group_name,
                                  const std::string& coef_name,
                                  QCDOrder order,
                                  ContributionType contribution,
                                  WilsonBasis basis = WilsonBasis::B_STANDARD)
{
    WGroupId group_id = GroupMapper::id_of(group_name);
    WCoefId coef_id = WCoefMapper::id_of(coef_name);

    return wilson.getFR(group_id, coef_id, order, contribution, basis);
}

} // namespace

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);

    // Make sure all builtin dynamic registries are populated before resolving strings.
    // HyperisoMaster usually does this during init, but this call is safe and idempotent.
    init_all_builtins();

    WilsonInterface wilson;

    // Build builtin Wilson groups from runtime strings.
    // These are the canonical names currently stored in GroupMapper.
    WilsonBuildConfig builtin_cfg(
        resolve_groups({"BCoefficients", "BPrimeCoefficients", "BScalarCoefficients"}),
        81.0,
        4.8,
        QCDOrder::LO
    );
    wilson.build(builtin_cfg);

    std::cout << "\nBuiltin coefficients queried through WGroupId/WCoefId overloads:\n";
    for (const std::string& coef_name : {"C7", "C9", "C10"}) {
        auto value = get_full_running_by_name(
            wilson,
            "BCoefficients",
            coef_name,
            QCDOrder::LO,
            ContributionType::TOTAL,
            WilsonBasis::B_STANDARD
        );
        std::cout << "BCoefficients::" << coef_name << " full running LO = " << value << "\n";
    }

    // True custom group + true custom coefficients.
    // The mapper creates runtime ids and the interface receives lambdas for the calculation.
    GroupMapper::register_custom("DEV_DYNAMIC_GROUP", {"dev-dynamic-group"});
    WCoefMapper::register_custom("C_DEV_A", {"cdevA"}, std::pair<int, int>{990100, 1});
    WCoefMapper::register_custom("C_DEV_B", {"cdevB"}, std::pair<int, int>{990100, 2});

    WGroupId custom_group = GroupMapper::id_of("dev-dynamic-group");
    WCoefId c_a = WCoefMapper::id_of("cdevA");
    WCoefId c_b = WCoefMapper::id_of("cdevB");

    CustomWilsonCoefficientConfig coef_a(c_a);
    coef_a.set_matching(
        QCDOrder::LO,
        {},
        [](const ParamSrc&) { return 1.50; },
        ContributionType::SM
    );

    CustomWilsonCoefficientConfig coef_b(c_b);
    coef_b.set_matching(
        QCDOrder::LO,
        {},
        [](const ParamSrc&) { return -0.25; },
        ContributionType::SM
    );

    CustomWilsonGroupConfig custom_cfg(custom_group);
    custom_cfg.matching_scale = 81.0;
    custom_cfg.hadronic_scale = 4.8;
    custom_cfg.order = QCDOrder::LO;
    custom_cfg.contribution = ContributionType::SM;
    custom_cfg.add_coefficient(coef_a).add_coefficient(coef_b);

    // Running lambda: here a toy identity-like evolution with a simple rescaling.
    // The input map is indexed by QCDOrder and WCoefId, so it supports custom ids directly.
    custom_cfg.set_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::LO,
        {},
        [c_a, c_b](const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
                   const BlockSrc&) {
            std::unordered_map<WCoefId, scalar_t> out;
            out[c_a] = 0.90 * matching.at(QCDOrder::LO).at(c_a);
            out[c_b] = 0.90 * matching.at(QCDOrder::LO).at(c_b);
            return out;
        }
    );

    wilson.addCustomWilsonGroup(custom_cfg);

    std::cout << "\nCustom coefficients queried through dynamic ids:\n";
    std::cout << "DEV_DYNAMIC_GROUP::C_DEV_A matching LO = "
              << wilson.getFM(custom_group, c_a, QCDOrder::LO, ContributionType::SM) << "\n";
    std::cout << "DEV_DYNAMIC_GROUP::C_DEV_A running LO = "
              << wilson.getFR(custom_group, c_a, QCDOrder::LO, ContributionType::SM) << "\n";
    std::cout << "DEV_DYNAMIC_GROUP::C_DEV_B running LO = "
              << wilson.getFR(custom_group, c_b, QCDOrder::LO, ContributionType::SM) << "\n";

    return 0;
}
