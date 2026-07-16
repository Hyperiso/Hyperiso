#include <cassert>
#include <iostream>

#include "GroupDefinition.h"
#include "CoefficientGroupBuilder.h"
#include "WilsonCoefficientRegistry.h"
#include "WilsonGroup.h"
#include "Wilson.h"

#include "stubs.hpp"

int main() {
    std::cout << "== GroupDefinition INTEGRATION ==\n";

    CoefficientRegistry reg;

    register_BScalar(reg);

    CoefficientGroupBuilder builder{reg};

    {
        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::BScalar)
        );

        const auto& def = GroupDefinitions::get(ctx.group_id);
        auto grp = builder.build(ctx);
        assert(grp && "BScalar(SM) group must be buildable");

        assert(grp->get_group_id() == def.id);
        assert(grp->get_type() == ContributionType::SM);

        auto bases = grp->get_bases();
        assert(bases.find(WilsonBasis::B_STANDARD) != bases.end());

        auto sources_LO = grp->get_sources(QCDOrder::LO, WilsonBasis::B_STANDARD);
        auto itW = sources_LO.find(ParameterType::WILSON);
        assert(itW != sources_LO.end());

        auto itBSM = sources_LO.find(ParameterType::BSM);
        assert(itBSM == sources_LO.end());
    }

    {
        auto ctx = make_ctx(
            Model::SUSY,
            Backend::Builtin,
            ContributionType::BSM,
            GroupMapper::to_id(WGroup::BScalar)
        );

        const auto& def = GroupDefinitions::get(ctx.group_id);
        auto grp = builder.build(ctx);
        assert(grp && "BScalar(SUSY) group must be buildable");

        assert(grp->get_group_id() == def.id);
        assert(grp->get_type() == ContributionType::BSM);

        auto bases = grp->get_bases();
        assert(bases.find(WilsonBasis::B_STANDARD) != bases.end());

        auto sources_LO = grp->get_sources(QCDOrder::LO, WilsonBasis::B_STANDARD);

        auto itW = sources_LO.find(ParameterType::WILSON);
        assert(itW != sources_LO.end());

        auto itSM = sources_LO.find(ParameterType::SM);
        assert(itSM != sources_LO.end());

        auto itBSM = sources_LO.find(ParameterType::BSM);
        assert(itBSM != sources_LO.end()
               && "SUSY hook for BScalar must add BSM sources");

        const auto& bsmBlocks = itBSM->second;
        assert(!bsmBlocks.empty());
    }

    std::cout << " GroupDefinition integration suite passed.\n";
    return 0;
}
