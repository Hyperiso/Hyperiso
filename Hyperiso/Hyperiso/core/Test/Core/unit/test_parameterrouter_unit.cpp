#include <cassert>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <algorithm>

#include "ParameterRouter.h"
#include "General.h"
#include "Utils.h"

int main() {
    std::cout << "== Running ParameterRouter unit tests ==\n";

    {
        std::vector<BlockName> in = {"MASS", "GAUGE", "FOO_UNKNOWN"};
        auto out = ParameterBlockRepartition::filter_custom_blocks(in);

        std::unordered_set<BlockName> s(out.begin(), out.end());
        assert(s.count("MASS") == 1);
        assert(s.count("GAUGE") == 1);
        assert(s.count("FOO_UNKNOWN") == 1);
    }


    {
        auto t1 = ParamRouter::GetType("MASS", LhaID(25));
        auto t2 = ParamRouter::GetType("MASS", LhaID(35));
        assert(t1 == ParameterType::SM);
        assert(t2 == ParameterType::BSM);

        auto g1 = ParamRouter::GetType("GAUGE", LhaID(2));
        auto g2 = ParamRouter::GetType("GAUGE", LhaID(99));
        assert(g1 == ParameterType::SM);
        assert(g2 == ParameterType::BSM);
    }


    {
        auto v_mass = ParamRouter::GetType("MASS");

        std::unordered_set<ParameterType> mass_types(v_mass.begin(), v_mass.end());
        assert(mass_types.count(ParameterType::SM) == 1);
        assert(mass_types.count(ParameterType::BSM) == 1);

        auto v_alpha = ParamRouter::GetType("ALPHA");

        std::unordered_set<ParameterType> alpha_types(v_alpha.begin(), v_alpha.end());
        assert(alpha_types.size() == 1 && alpha_types.count(ParameterType::BSM) == 1);
    }


    {
        auto sm_blocks = ParamRouter::GetOwnedBlocks(ParameterType::SM);
        assert(sm_blocks.count("MASS") == 1);
        assert(sm_blocks.count("GAUGE") == 1);
        assert(sm_blocks.count("SMINPUTS") == 1);

        auto bsm_blocks = ParamRouter::GetOwnedBlocks(ParameterType::BSM);
        assert(bsm_blocks.count("MASS") == 1);
        assert(bsm_blocks.count("ALPHA") == 1);
    }

    std::cout << "\n ParameterRouter unit suite passed.\n";
    return 0;
}
