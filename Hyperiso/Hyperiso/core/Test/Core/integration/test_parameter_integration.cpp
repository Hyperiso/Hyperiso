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

int main() {
    std::cout << "== Running Parameter integration tests ==\n";

    auto pA = make_param(1.0, "A", 1);
    auto pB = make_param(2.0, "B", 2);

    auto depSum = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "SUM", 10},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {pA->get_id(), pA},
            {pB->get_id(), pB}
        },
        [](const auto& srcs, std::shared_ptr<DependentParameter> self) {
            double s = 0.0; for (auto& [_, p] : srcs.raw()) s += p->get_val();
            self->set_expected(s);
        }
    );
    depSum->init();

    auto depScaled = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "SCALED", 20},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {depSum->get_id(), depSum}
        },
        [](const auto& srcs, std::shared_ptr<DependentParameter> self) {
            double base = srcs.raw().begin()->second->get_val();
            self->set_expected(3.0 * base);
        }
    );
    depScaled->init();

    depSum->update();
    depScaled->update();
    assert(std::abs(depSum->get_val()    - (1.0 + 2.0)) < 1e-12);
    assert(std::abs(depScaled->get_val() - 3.0 * (1.0 + 2.0)) < 1e-12);

    pB->set_expected(5.0);
    assert(std::abs(depSum->get_val()    - (1.0 + 5.0)) < 1e-12);
    assert(std::abs(depScaled->get_val() - 3.0 * (1.0 + 5.0)) < 1e-12);

    depSum->freeze();
    pA->set_expected(7.0);

    assert(std::abs(depSum->get_val()    - (1.0 + 5.0)) < 1e-12);
    assert(std::abs(depScaled->get_val() - 3.0 * (1.0 + 5.0)) < 1e-12);

    depSum->unfreeze();
    assert(std::abs(depSum->get_val()    - (7.0 + 5.0)) < 1e-12);
    assert(std::abs(depScaled->get_val() - 3.0 * (7.0 + 5.0)) < 1e-12);

    depScaled->clear_above();
    pB->set_expected(6.0);

    assert(std::abs(depSum->get_val()    - (7.0 + 6.0)) < 1e-12);
    assert(std::abs(depScaled->get_val() - 3.0 * (7.0 + 5.0)) < 1e-12);

    std::cout << "\n Parameter integration suite passed.\n";
    return 0;
}
