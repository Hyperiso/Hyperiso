#include <cassert>
#include <iostream>

#include "GroupDefinition.h"
#include "WilsonGroup.h"
#include "Wilson.h"

int main() {
    std::cout << "== GroupDefinition UNIT ==\n";


    {
        const auto& defB = GroupDefinitions::B();

        auto expected_id = GroupMapper::to_id(WGroup::B);
        assert(defB.id == expected_id);

        assert(!defB.members.empty());

        auto it_basis_std = defB.sources.find(WilsonBasis::B_STANDARD);
        assert(it_basis_std != defB.sources.end());

        const auto& per_order_std = it_basis_std->second;
        assert(per_order_std.find(QCDOrder::LO)   != per_order_std.end());
        assert(per_order_std.find(QCDOrder::NLO)  != per_order_std.end());
        assert(per_order_std.find(QCDOrder::NNLO) != per_order_std.end());

        const auto& lo_cgs_std = per_order_std.at(QCDOrder::LO);
        auto it_wilson = lo_cgs_std.sources.find(ParameterType::WILSON);
        assert(it_wilson != lo_cgs_std.sources.end());

        const auto& namesW = it_wilson->second;
        bool has_placeholder = false;
        for (const auto& n : namesW) {
            if (n == MATCHING_BLOCK_PLACEHOLDER) {
                has_placeholder = true;
                break;
            }
        }
        assert(has_placeholder && "B group LO WILSON sources must contain MATCHING_BLOCK_PLACEHOLDER");

        auto it_basis_trad = defB.sources.find(WilsonBasis::B_TRADITIONAL);
        assert(it_basis_trad != defB.sources.end());
        const auto& per_order_trad = it_basis_trad->second;
        assert(per_order_trad.find(QCDOrder::LO) != per_order_trad.end());
    }

    {
        const auto& defBp = GroupDefinitions::BPrime();
        auto expected_id = GroupMapper::to_id(WGroup::BPrime);
        assert(defBp.id == expected_id);
        assert(!defBp.members.empty());

        auto it_basis = defBp.sources.find(WilsonBasis::B_STANDARD);
        assert(it_basis != defBp.sources.end());
        const auto& per_order = it_basis->second;
        assert(per_order.find(QCDOrder::LO) != per_order.end());

        const auto& lo_cgs = per_order.at(QCDOrder::LO);
        auto it_wilson = lo_cgs.sources.find(ParameterType::WILSON);
        assert(it_wilson != lo_cgs.sources.end());
        const auto& namesW = it_wilson->second;
        bool has_placeholder = false;
        for (const auto& n : namesW) {
            if (n == MATCHING_BLOCK_PLACEHOLDER) {
                has_placeholder = true;
                break;
            }
        }
        assert(has_placeholder && "BPrime LO WILSON sources must contain MATCHING_BLOCK_PLACEHOLDER");
    }

    {
        const auto& defBs = GroupDefinitions::BScalar();
        auto expected_id = GroupMapper::to_id(WGroup::BScalar);
        assert(defBs.id == expected_id);
        assert(!defBs.members.empty());

        auto it_basis = defBs.sources.find(WilsonBasis::B_STANDARD);
        assert(it_basis != defBs.sources.end());
        const auto& per_order = it_basis->second;
        assert(per_order.find(QCDOrder::LO) != per_order.end());
        assert(per_order.find(QCDOrder::NLO) != per_order.end());

        auto it_susy = defBs.setup.find(Model::SUSY);
        assert(it_susy != defBs.setup.end());
        assert(!it_susy->second.empty()
               && "BScalar should register a SUSY setup hook");
    }

    {
        const auto& defK = GroupDefinitions::K();
        auto expected_id = GroupMapper::to_id(WGroup::K);
        assert(defK.id == expected_id);
        assert(!defK.members.empty());

        auto it_basis = defK.sources.find(WilsonBasis::B_STANDARD);
        assert(it_basis != defK.sources.end());
        const auto& per_order = it_basis->second;
        assert(per_order.find(QCDOrder::LO) != per_order.end());

        const auto& lo_cgs = per_order.at(QCDOrder::LO);
        auto it_wilson = lo_cgs.sources.find(ParameterType::WILSON);
        assert(it_wilson != lo_cgs.sources.end());
        const auto& namesW = it_wilson->second;
        bool has_placeholder = false;
        for (const auto& n : namesW) {
            if (n == MATCHING_BLOCK_PLACEHOLDER) {
                has_placeholder = true;
                break;
            }
        }
        assert(has_placeholder && "K LO WILSON sources must contain MATCHING_BLOCK_PLACEHOLDER");
    }

    {

        WGroupId custom_id{};

        GroupDefinition def;
        def.id = custom_id;
        def.members.clear();
        def.sources.clear();

        GroupDefinitions::register_custom(def);

        assert(GroupDefinitions::has_custom(custom_id));

        const auto& got = GroupDefinitions::get_custom(custom_id);
        assert(got.id == custom_id);

        const auto& got2 = GroupDefinitions::get(custom_id);
        assert(&got == &got2);
    }

    std::cout << " GroupDefinition unit suite passed.\n";
    return 0;
}
