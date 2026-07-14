#include <cassert>
#include <iostream>

#include "WilsonCoefficientRegistry.h"
#include "GroupDefinition.h" 
#include "Wilson.h"  
#include "BWilson.h"
#include "KWilson.h"  
#include "MartyWilson.h"
#include "BWilsonSUSY.h"

#include "stubs.hpp"

int main() {
    std::cout << "== WilsonCoefficientRegistry UNIT ==\n";

    {
        CoefficientRegistry empty;
        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::B)
        );

        bool threw = false;
        try {
            (void)empty.create(ctx, WCoef::C1);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw && "Registry must throw when no factory is registered for a coefficient");
    }

    {
        CoefficientRegistry reg;
        register_B(reg);

        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::B)
        );

        auto c1 = reg.create(ctx, WCoef::C1);
        assert(c1 && "C1(SM, Builtin) should be creatable");
        auto asC1 = dynamic_cast<C1*>(c1.get());
        assert(asC1 && "C1(SM, Builtin) should be instance of C1");
    }

    {
        CoefficientRegistry reg;
        register_B(reg);

        auto ctx = make_ctx(
            Model::SUSY,
            Backend::Builtin,
            ContributionType::BSM,
            GroupMapper::to_id(WGroup::B)
        );

        auto c1_susy = reg.create(ctx, WCoef::C1);
        assert(c1_susy && "C1(SUSY, Builtin) should be creatable");
        auto asC1_susy = dynamic_cast<C1_susy*>(c1_susy.get());
        assert(asC1_susy && "C1(SUSY, Builtin) should be instance of C1_susy");
    }

    {
        CoefficientRegistry reg;
        register_B(reg);

        auto ctx = make_ctx(
            Model::SUSY,
            Backend::Marty, 
            ContributionType::BSM,
            GroupMapper::to_id(WGroup::B)
        );

        auto c1_marty = reg.create(ctx, WCoef::C1);
        assert(c1_marty && "C1(SUSY, Marty, BSM) should be creatable through MARTY");
        auto as_marty = dynamic_cast<MartyWilson*>(c1_marty.get());
        assert(as_marty && "A BSM MARTY request should use the MARTY coefficient factory");
    }

    {
        CoefficientRegistry reg;
        register_K(reg);

        auto ctx = make_ctx(
            Model::SUSY,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::K)
        );

        auto ck9 = reg.create(ctx, WCoef::CK9);
        assert(ck9 && "CK9(SUSY, Builtin) should fall back to CK9(SM, Builtin)");
        auto asCK9 = dynamic_cast<CK9*>(ck9.get());
        assert(asCK9 && "Fallback model->SM should return CK9 (SM implementation)");
    }

    std::cout << " WilsonCoefficientRegistry unit suite passed.\n";
    return 0;
}
