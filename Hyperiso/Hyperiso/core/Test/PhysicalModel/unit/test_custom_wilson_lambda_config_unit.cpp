#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <utility>

#include "CustomWilsonLambda.h"
#include "Parameter.h"
#include "SourcesView.h"
#include "mapper_hub.hpp"

static bool approx(double a, double b, double eps = 1e-12) {
    return std::fabs(a - b) <= eps;
}

int main() {
    std::cout << "== CustomWilson lambda config UNIT ==\n";

    init_all_builtins();

    GroupMapper::register_custom("UNIT_LAMBDA_WILSON_GROUP", {"unit-lambda-wilson"});
    WCoefMapper::register_custom("C_UNIT_LAMBDA", {"c-unit-lambda"}, std::pair<int, int>{920001, 1});

    WGroupId group = GroupMapper::id_of("unit-lambda-wilson");
    WCoefId coeff = WCoefMapper::id_of("c-unit-lambda");

    ParamId mt(ParameterType::SM, "SMINPUTS", LhaID(6));
    ParamId fb(ParameterType::FLAVOR, "FCONST", LhaID(531, 1));

    CustomWilsonCoefficientConfig coef_cfg(coeff);
    coef_cfg.set_matching(
        QCDOrder::LO,
        {mt, fb},
        [mt, fb](const ParamSrc& src) {
            return 0.1 * src.get_val(mt) + 2.0 * src.get_val(fb);
        },
        ContributionType::SM
    );

    assert(coef_cfg.id == coeff);
    assert(coef_cfg.matching.count(QCDOrder::LO) == 1);
    assert(coef_cfg.matching.at(QCDOrder::LO).sources.count(mt) == 1);
    assert(coef_cfg.matching.at(QCDOrder::LO).sources.count(fb) == 1);

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> params;
    params.emplace(mt, std::make_shared<Parameter>(mt, 173.0, 0.0, 0.0));
    params.emplace(fb, std::make_shared<Parameter>(fb, 0.230, 0.0, 0.0));
    ParamSrc src(params, "unit custom Wilson sources");

    const double matching_value = std::real(coef_cfg.matching.at(QCDOrder::LO).compute(src));
    assert(approx(matching_value, 0.1 * 173.0 + 2.0 * 0.230));

    CustomWilsonGroupConfig group_cfg(group);
    group_cfg.matching_scale = 81.0;
    group_cfg.hadronic_scale = 4.8;
    group_cfg.order = QCDOrder::LO;
    group_cfg.add_coefficient(coef_cfg);

    group_cfg.set_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::LO,
        {},
        [coeff](const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
                const BlockSrc&) {
            std::unordered_map<WCoefId, scalar_t> out;
            out[coeff] = 1.5 * matching.at(QCDOrder::LO).at(coeff);
            return out;
        }
    );

    assert(group_cfg.group == group);
    assert(group_cfg.coefficients.size() == 1);
    assert(group_cfg.running.count(WilsonBasis::B_STANDARD) == 1);
    assert(group_cfg.running.at(WilsonBasis::B_STANDARD).count(QCDOrder::LO) == 1);

    std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> matching;
    matching[QCDOrder::LO][coeff] = scalar_t{2.0, 0.0};
    std::unordered_map<std::string, std::shared_ptr<Block>> blocks;
    BlockSrc block_src(blocks, "unit custom Wilson blocks");

    auto running = group_cfg.running.at(WilsonBasis::B_STANDARD).at(QCDOrder::LO).func(matching, block_src);
    assert(running.count(coeff) == 1);
    assert(approx(std::real(running.at(coeff)), 3.0));

    auto identity = custom_identity_running(matching, block_src);
    assert(identity.count(coeff) == 1);
    assert(approx(std::real(identity.at(coeff)), 2.0));

    std::cout << "UNIT OK\n";
    return 0;
}
