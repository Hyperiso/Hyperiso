#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Parameter.h"
#include "DependentParameter.h"
#include "Include.h"
#include "SourcesView.h"

static std::shared_ptr<Parameter> make_param(double v, const char* block, int code) {
    auto p = std::make_shared<Parameter>(ParamId{ParameterType::SM, block, code}, v, 0.0, 0.0);
    p->set_expected(v);
    return p;
}

static void test_parameter_chain_detach_reattach() {
    std::cout << "\n-- [Integration] Parameter chain detach/reattach --\n";

    auto pA = make_param(1.0, "A2", 1);
    auto pB = make_param(2.0, "B2", 2);

    int runsSum = 0;
    auto depSum = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "SUM2", 10},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {pA->get_id(), pA},
            {pB->get_id(), pB}
        },
        [&](const auto& srcs, std::shared_ptr<DependentParameter> self) {
            ++runsSum;
            double s = 0.0;
            for (auto& [_, p] : srcs.raw()) s += p->get_val();
            self->set_expected_silent(s);
        }
    );
    depSum->init();

    int runsScaled = 0;
    auto depScaled = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "SCALED2", 20},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {depSum->get_id(), depSum}
        },
        [&](const auto& srcs, std::shared_ptr<DependentParameter> self) {
            ++runsScaled;
            double base = srcs.raw().begin()->second->get_val();
            self->set_expected_silent(3.0 * base);
        }
    );
    depScaled->init();

    assert(std::abs(depScaled->get_val() - 9.0) < 1e-12);
    assert(runsSum == 1);
    assert(runsScaled == 1);

    depSum->detach();
    pA->set_expected(10.0);

    assert(std::abs(depSum->get_val() - 3.0) < 1e-12);
    assert(std::abs(depScaled->get_val() - 9.0) < 1e-12);
    assert(runsSum == 1);
    assert(runsScaled == 1);

    depSum->reattach();

    assert(std::abs(depSum->get_val() - 12.0) < 1e-12);
    assert(std::abs(depScaled->get_val() - 36.0) < 1e-12);
    assert(runsSum == 2);
    assert(runsScaled == 2);

    depScaled->detach();
    pB->set_expected(5.0);

    assert(std::abs(depSum->get_val() - 15.0) < 1e-12);
    assert(std::abs(depScaled->get_val() - 36.0) < 1e-12);
    assert(runsSum == 3);
    assert(runsScaled == 2);

    depScaled->reattach();
    assert(std::abs(depScaled->get_val() - 45.0) < 1e-12);
    assert(runsScaled == 3);
}

int main() {
    std::cout << "== Running Parameter integration tests ==\n";

    auto pA = make_param(1.0, "A", 1);
    auto pB = make_param(2.0, "B", 2);

    int runsSum = 0;
    auto depSum = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "SUM", 10},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {pA->get_id(), pA},
            {pB->get_id(), pB}
        },
        [&](const auto& srcs, std::shared_ptr<DependentParameter> self) {
            ++runsSum;
            double s = 0.0;
            for (auto& [_, p] : srcs.raw()) s += p->get_val();
            self->set_expected(s);
        }
    );
    depSum->init();

    int runsScaled = 0;
    auto depScaled = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "SCALED", 20},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {depSum->get_id(), depSum}
        },
        [&](const auto& srcs, std::shared_ptr<DependentParameter> self) {
            ++runsScaled;
            double base = srcs.raw().begin()->second->get_val();
            self->set_expected(3.0 * base);
        }
    );
    depScaled->init();

    assert(runsSum == 0);
    assert(runsScaled == 0);

    assert(std::abs(depScaled->get_val() - 3.0 * (1.0 + 2.0)) < 1e-12);
    assert(std::abs(depSum->get_val() - (1.0 + 2.0)) < 1e-12);
    assert(runsSum == 1);
    assert(runsScaled == 1);

    pB->set_expected(5.0);
    assert(runsSum == 1);
    assert(runsScaled == 1);

    assert(std::abs(depScaled->get_val() - 3.0 * (1.0 + 5.0)) < 1e-12);
    assert(std::abs(depSum->get_val() - (1.0 + 5.0)) < 1e-12);
    assert(runsSum == 2);
    assert(runsScaled == 2);

    depSum->freeze();
    pA->set_expected(7.0);

    assert(std::abs(depSum->get_val() - (1.0 + 5.0)) < 1e-12);

    assert(std::abs(depScaled->get_val() - 3.0 * (1.0 + 5.0)) < 1e-12);

    depSum->unfreeze();
    assert(std::abs(depScaled->get_val() - 3.0 * (7.0 + 5.0)) < 1e-12);
    assert(std::abs(depSum->get_val() - (7.0 + 5.0)) < 1e-12);

    depScaled->clear_above();

    pB->set_expected(6.0);

    assert(std::abs(depSum->get_val() - (7.0 + 6.0)) < 1e-12);
    assert(std::abs(depScaled->get_val() - 3.0 * (7.0 + 5.0)) < 1e-12);

    test_parameter_chain_detach_reattach();
    
    std::cout << "\n Parameter integration suite passed.\n";
    return 0;
}
