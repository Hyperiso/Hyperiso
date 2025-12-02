#include <cassert>
#include <iostream>

#include "WilsonCoefficientRegistry.h"
#include "GroupDefinition.h"
#include "Wilson.h"

#include "stubs.hpp"

int main() {
    std::cout << "== WilsonCoefficientRegistry INTEGRATION ==\n";

    CoefficientRegistry reg;

    register_B(reg);
    register_BPrime(reg);
    register_BScalar(reg);
    register_CC_bc(reg);
    register_CC_bu(reg);
    register_CC_cs(reg);
    register_CC_cd(reg);
    register_CC_su(reg);
    register_CC_du(reg);
    register_MesonMixing(reg);
    register_K(reg);

    {
        auto ctx_B = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::B)
        );
        auto c7 = reg.create(ctx_B, WCoef::C7);
        assert(c7 && "C7(SM, Builtin) should be creatable");

        auto ctx_Bp = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::BPrime)
        );
        auto cp1 = reg.create(ctx_Bp, WCoef::CP1);
        assert(cp1 && "CP1(SM, Builtin) should be creatable");

        auto ctx_Bs = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::BScalar)
        );
        auto cq1 = reg.create(ctx_Bs, WCoef::CQ1);
        assert(cq1 && "CQ1(SM, Builtin) should be creatable");

        auto ctx_cc_bc = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::CC_bc)
        );
        auto cv1_bc = reg.create(ctx_cc_bc, WCoef::C_V1_bc);
        assert(cv1_bc && "C_V1_bc(SM, Builtin) should be creatable");

        auto ctx_K = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::K)
        );
        auto ck9 = reg.create(ctx_K, WCoef::CK9);
        assert(ck9 && "CK9(SM, Builtin) should be creatable");
    }

    {
        auto ctx_B_susy = make_ctx(
            Model::SUSY,
            Backend::Builtin,
            ContributionType::BSM,
            GroupMapper::to_id(WGroup::B)
        );
        auto c7_susy = reg.create(ctx_B_susy, WCoef::C7);
        assert(c7_susy && "C7(SUSY, Builtin) should be creatable");

        auto ctx_B_thdm = make_ctx(
            Model::THDM,
            Backend::Builtin,
            ContributionType::BSM,
            GroupMapper::to_id(WGroup::B)
        );
        auto c7_thdm = reg.create(ctx_B_thdm, WCoef::C7);
        assert(c7_thdm && "C7(THDM, Builtin) should be creatable");
    }

    {
        auto ctx_K_susy = make_ctx(
            Model::SUSY,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::K)
        );
        auto ck9_susy = reg.create(ctx_K_susy, WCoef::CK9);
        assert(ck9_susy && "CK9(SUSY, Builtin) should fall back to CK9(SM, Builtin)");

        auto ctx_K_thdm = make_ctx(
            Model::THDM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::K)
        );
        auto ck9_thdm = reg.create(ctx_K_thdm, WCoef::CK9);
        assert(ck9_thdm && "CK9(THDM, Builtin) should fall back to CK9(SM, Builtin)");
    }

    std::cout << " WilsonCoefficientRegistry integration suite passed.\n";
    return 0;
}
