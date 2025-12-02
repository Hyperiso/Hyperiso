#include <cassert>
#include <iostream>

#include "CoefficientGroupBuilder.h"
#include "WilsonCoefficientRegistry.h"
#include "GroupDefinition.h"
#include "WilsonGroup.h"
#include "Wilson.h"
#include "GroupMapper.h"
#include "wcoef_ids.hpp"

#include "stubs.hpp"

int main() {
    std::cout << "== CoefficientGroupBuilder UNIT ==\n";

    CoefficientRegistry reg;
    register_B(reg);

    CoefficientGroupBuilder builder{reg};

    auto ctx = make_ctx(
        Model::SM,
        Backend::Builtin,
        ContributionType::SM,
        GroupMapper::to_id(WGroup::B)
    );

    const auto& def = GroupDefinitions::get(ctx.group_id);

    auto grp = builder.build(ctx);
    assert(grp && "Builder must return a non-null CoefficientGroup");

    assert(grp->get_group_id() == def.id);

    assert(grp->get_type() == ContributionType::SM);

    {
        for (auto c : def.members) {
            auto name = WCoefMapper::str(c);
            auto it   = grp->find(name);
            assert(it != grp->end() && "All group members must be present in the built CoefficientGroup");
            assert(it->second && "Coefficient pointers must be non-null");
        }
    }

    {
        auto bases = grp->get_bases();
        assert(bases.find(WilsonBasis::B_STANDARD) != bases.end()
               && "B_STANDARD basis should be present for B group");

        auto sources_LO = grp->get_sources(QCDOrder::LO, WilsonBasis::B_STANDARD);
        auto itW = sources_LO.find(ParameterType::WILSON);
        assert(itW != sources_LO.end() && "WILSON sources must be present for B group / B_STANDARD / LO");

        const auto& names = itW->second;

        std::string matching_block = GroupMapper::str(def.id, ScaleType::MATCHING);

        bool has_matching_block = false;
        bool has_placeholder = false;
        for (const auto& n : names) {
            if (n == matching_block) has_matching_block = true;
            if (n == MATCHING_BLOCK_PLACEHOLDER) has_placeholder = true;
        }

        assert(has_matching_block && "Matching block name must appear in WILSON sources after build()");
        assert(!has_placeholder && "MATCHING_BLOCK_PLACEHOLDER must have been replaced");
    }

    std::cout << " CoefficientGroupBuilder unit suite passed.\n";
    return 0;
}
