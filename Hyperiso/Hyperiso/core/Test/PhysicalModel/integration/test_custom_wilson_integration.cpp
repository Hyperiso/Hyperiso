#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CustomWilsonLambda.h"
#include "IBlockComposer.h"
#include "registry_init.hpp"
#include "SourcesView.h"
#include "Parameter.h"
#include "Block.h"

#include "stubs.hpp"

namespace {

struct RecordingComposer : IBlockComposer {
    struct ParamCall {
        ParamId target;
        std::unordered_set<ParamId> sources;
        DepParamUpdateFunc update;
    };

    std::vector<ParamCall> param_calls;

    void compose_block(const std::string&, const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& id, const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc& update) override {
        param_calls.push_back({id, src, update});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

std::shared_ptr<Parameter> param(const ParamId& id, double value) {
    return std::make_shared<Parameter>(id, value, 0.0, 0.0);
}

} // namespace

int main() {
    std::cout << "== Custom Wilson lambda INTEGRATION ==\n";

    init_all_builtins();

    GroupMapper::register_custom("TEST_INTEGRATION_WILSON_GROUP", {"test_integration_group"});
    WCoefMapper::register_custom("C_TEST_INT_A", {"c_test_int_a"}, std::make_pair(930201, 1));
    WCoefMapper::register_custom("C_TEST_INT_B", {"c_test_int_b"}, std::make_pair(930201, 2));

    const WGroupId group = GroupMapper::id_of("test_integration_group");
    const WCoefId c_a = WCoefMapper::id_of("c_test_int_a");
    const WCoefId c_b = WCoefMapper::id_of("c_test_int_b");

    const ParamId p_a{ParameterType::SM, "SMINPUTS", LhaID(6)};
    const ParamId p_b{ParameterType::FLAVOR, "FCONST", LhaID(531, 1)};

    CustomWilsonCoefficientConfig coeff_a(c_a);
    coeff_a.set_matching(QCDOrder::LO, {p_a}, [p_a](const ParamSrc& src) {
        return src.get_val(p_a) / 100.0;
    }, ContributionType::BSM);

    CustomWilsonCoefficientConfig coeff_b(c_b);
    coeff_b.set_matching(QCDOrder::LO, {p_b}, [p_b](const ParamSrc& src) {
        return 2.0 * src.get_val(p_b);
    }, ContributionType::BSM);

    CustomWilsonGroupConfig cfg(group);
    cfg.display_name = "TEST_INTEGRATION_WILSON_GROUP";
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
            out[c_b] = matching.at(QCDOrder::LO).at(c_a) + matching.at(QCDOrder::LO).at(c_b);
            return out;
        }
    );

    auto ctx = make_ctx(Model::SM, Backend::Builtin, ContributionType::BSM, group);
    auto recorder = std::make_shared<RecordingComposer>();
    ctx.adapters.iblock_c = recorder;

    auto custom_group = make_custom_wilson_group(cfg, ctx.adapters);
    custom_group->finalize(QCDOrder::LO);

    assert(custom_group->get_order() == QCDOrder::LO);
    assert(recorder->param_calls.size() == 2);

    bool saw_a = false;
    bool saw_b = false;
    for (const auto& call : recorder->param_calls) {
        if (call.target == ParamId{WilsonBlockNames::matching(group), WCoefMapper::flha_full(c_a, QCDOrder::LO, ContributionType::BSM)}) {
            assert(call.sources.count(p_a) == 1);
            saw_a = true;
        }
        if (call.target == ParamId{WilsonBlockNames::matching(group), WCoefMapper::flha_full(c_b, QCDOrder::LO, ContributionType::BSM)}) {
            assert(call.sources.count(p_b) == 1);
            saw_b = true;
        }
    }
    assert(saw_a);
    assert(saw_b);

    // Vérifie que les lambdas réellement composées par le groupe calculent les valeurs attendues.
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> values{
        {p_a, param(p_a, 173.0)},
        {p_b, param(p_b, 0.23)}
    };

    auto f_a = custom_group->at(WCoefMapper::str(c_a))->get_func(QCDOrder::LO);
    auto f_b = custom_group->at(WCoefMapper::str(c_b))->get_func(QCDOrder::LO);
    assert(std::abs(std::real(f_a(ParamSrc(values))) - 1.73) < 1e-12);
    assert(std::abs(std::real(f_b(ParamSrc(values))) - 0.46) < 1e-12);

    std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> matching;
    matching[QCDOrder::LO][c_a] = f_a(ParamSrc(values));
    matching[QCDOrder::LO][c_b] = f_b(ParamSrc(values));
    std::unordered_map<std::string, std::shared_ptr<Block>> blocks;
    auto running = custom_group->get_func(QCDOrder::LO, WilsonBasis::B_STANDARD)(matching, BlockSrc(blocks));
    assert(std::abs(std::real(running.at(c_a)) - 1.73) < 1e-12);
    assert(std::abs(std::real(running.at(c_b)) - 2.19) < 1e-12);

    std::cout << " Custom Wilson lambda integration suite passed.\n";
    return 0;
}
