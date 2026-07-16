#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "CustomWilsonLambda.h"
#include "registry_init.hpp"
#include "SourcesView.h"
#include "Parameter.h"
#include "Block.h"

#include "stubs.hpp"

namespace {

std::shared_ptr<Parameter> param(const ParamId& id, double value) {
    return std::make_shared<Parameter>(id, value, 0.0, 0.0);
}

bool contains_id(const std::vector<WCoefId>& ids, const WCoefId& id) {
    return std::find(ids.begin(), ids.end(), id) != ids.end();
}

} // namespace

int main() {
    std::cout << "== Custom Wilson lambda UNIT ==\n";

    init_all_builtins();

    GroupMapper::register_custom("TEST_UNIT_WILSON_GROUP", {"test_unit_group"});
    WCoefMapper::register_custom("C_TEST_UNIT_A", {"c_test_unit_a"}, std::make_pair(930101, 1));
    WCoefMapper::register_custom("C_TEST_UNIT_B", {"c_test_unit_b"}, std::make_pair(930101, 2));

    const WGroupId group = GroupMapper::id_of("test_unit_group");
    const WCoefId c_a = WCoefMapper::id_of("c_test_unit_a");
    const WCoefId c_b = WCoefMapper::id_of("c_test_unit_b");

    const ParamId mtop{ParameterType::SM, "SMINPUTS", LhaID(6)};
    const ParamId alpha{ParameterType::SM, "SMINPUTS", LhaID(3)};

    CustomWilsonCoefficientConfig coeff_a(c_a);
    coeff_a.set_matching(
        QCDOrder::LO,
        {mtop},
        [mtop](const ParamSrc& src) { return src.get_val(mtop) / 100.0; },
        ContributionType::BSM
    );

    CustomWilsonCoefficientConfig coeff_b(c_b);
    coeff_b.set_matching(
        QCDOrder::LO,
        {alpha},
        [alpha](const ParamSrc& src) { return 10.0 * src.get_val(alpha); },
        ContributionType::BSM
    );

    // Teste directement le coefficient dynamique : id runtime, LHA complet et lambda utilisateur.
    CustomWilson direct(c_a, WilsonBlockNames::matching(group), ContributionType::BSM);
    direct.set_order_info(
        QCDOrder::LO,
        coeff_a.matching.at(QCDOrder::LO).sources,
        coeff_a.matching.at(QCDOrder::LO).compute,
        WCoefMapper::flha_full(c_a, QCDOrder::LO, ContributionType::BSM)
    );

    assert(direct.get_name() == WCoefMapper::str(c_a));
    assert(direct.get_storage_block() == WilsonBlockNames::matching(group));
    assert(direct.get_lhaid(QCDOrder::LO).to_string() ==
           WCoefMapper::flha_full(c_a, QCDOrder::LO, ContributionType::BSM).to_string());
    assert(direct.get_sources(QCDOrder::LO).count(mtop) == 1);

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> params{{mtop, param(mtop, 172.0)}};
    assert(std::abs(std::real(direct.get_func(QCDOrder::LO)(ParamSrc(params))) - 1.72) < 1e-12);

    CustomWilsonGroupConfig cfg(group);
    cfg.display_name = "TEST_UNIT_WILSON_GROUP";
    cfg.contribution = ContributionType::BSM;
    cfg.order = QCDOrder::LO;
    cfg.add_coefficient(coeff_a);
    cfg.add_coefficient(coeff_b);
    cfg.set_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::LO,
        {{ParameterType::WILSON, {WilsonBlockNames::matching(group)}}},
        [c_a, c_b](const auto& matching, const BlockSrc&) {
            std::unordered_map<WCoefId, scalar_t> out;
            out[c_a] = matching.at(QCDOrder::LO).at(c_a);
            out[c_b] = matching.at(QCDOrder::LO).at(c_b) + matching.at(QCDOrder::LO).at(c_a);
            return out;
        }
    );

    auto ctx = make_ctx(Model::SM, Backend::Builtin, ContributionType::BSM, group);
    auto custom_group = make_custom_wilson_group(cfg, ctx.adapters);

    assert(custom_group->get_group_id() == group);
    assert(custom_group->get_type() == ContributionType::BSM);
    assert(custom_group->size() == 2);
    assert(custom_group->find(WCoefMapper::str(c_a)) != custom_group->end());
    assert(custom_group->find(WCoefMapper::str(c_b)) != custom_group->end());
    assert(contains_id(custom_group->get_member_ids(), c_a));
    assert(contains_id(custom_group->get_member_ids(), c_b));

    auto sources = custom_group->get_sources(QCDOrder::LO, WilsonBasis::B_STANDARD);
    assert(sources.count(ParameterType::WILSON) == 1);
    assert(sources.at(ParameterType::WILSON).front() == WilsonBlockNames::matching(group));

    std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> matching;
    matching[QCDOrder::LO][c_a] = 1.25;
    matching[QCDOrder::LO][c_b] = 0.50;
    std::unordered_map<std::string, std::shared_ptr<Block>> blocks;
    auto running = custom_group->get_func(QCDOrder::LO, WilsonBasis::B_STANDARD)(matching, BlockSrc(blocks));
    assert(std::abs(std::real(running.at(c_a)) - 1.25) < 1e-12);
    assert(std::abs(std::real(running.at(c_b)) - 1.75) < 1e-12);

    std::cout << " Custom Wilson lambda unit suite passed.\n";
    return 0;
}
