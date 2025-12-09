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
    std::cout << "== CoefficientGroupBuilder INTEGRATION ==\n";

    CoefficientRegistry reg;

    register_B(reg);
    register_BPrime(reg);
    register_BScalar(reg);
    register_CC_bc(reg);
    register_K(reg);

    CoefficientGroupBuilder builder{reg};

    {
        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::B)
        );

        const auto& def = GroupDefinitions::get(ctx.group_id);
        auto grp = builder.build(ctx);
        assert(grp && "B group CoefficientGroup should be constructible");

        assert(grp->get_group_id() == def.id);
        assert(grp->get_type() == ContributionType::SM);

        for (auto c : def.members) {
            auto name = WCoefMapper::str(c);
            assert(grp->find(name) != grp->end());
        }
    }


    {
        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::BPrime)
        );

        const auto& def = GroupDefinitions::get(ctx.group_id);
        auto grp = builder.build(ctx);
        assert(grp && "BPrime group CoefficientGroup should be constructible");

        assert(grp->get_group_id() == def.id);

        for (auto c : def.members) {
            auto name = WCoefMapper::str(c);
            assert(grp->find(name) != grp->end());
        }
    }


    {
        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::CC_bc)
        );

        const auto& def = GroupDefinitions::get(ctx.group_id);
        auto grp = builder.build(ctx);
        assert(grp && "CC_bc group CoefficientGroup should be constructible");

        assert(grp->get_group_id() == def.id);

        for (auto c : def.members) {
            auto name = WCoefMapper::str(c);
            assert(grp->find(name) != grp->end());
        }
    }

    {
        auto ctx = make_ctx(
            Model::SUSY,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::K)
        );

        const auto& def = GroupDefinitions::get(ctx.group_id);
        auto grp = builder.build(ctx);
        assert(grp && "K group CoefficientGroup should be constructible even for non-SM models via fallback");

        assert(grp->get_group_id() == def.id);

        for (auto c : def.members) {
            auto name = WCoefMapper::str(c);
            assert(grp->find(name) != grp->end());
        }
    }

    std::cout << " CoefficientGroupBuilder integration suite passed.\n";
    return 0;
}
