#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Parameter.h"
#include "DependentParameter.h"
#include "SourcesView.h"
#include "Include.h" 

struct CountingParam : Parameter {
    int updates = 0;
    CountingParam() : Parameter(ParamId{ParameterType::SM, "TEST", 0}, 0.0, 0.0, 0.0) {}
    void update() override { ++updates; }
};

static void test_construct_and_std_and_mode() {
    std::cout << "\n-- [Unit] construct/std/mode --\n";

    ParamId pid{ParameterType::SM, "GAUGE", 3};
    Parameter p(pid, 0.1, 0.01, 0.02);

    assert(std::abs(p.get_val() - 0.1) < 1e-12);

    auto [st, sy] = p.get_std();
    assert(std::abs(st - 0.01) < 1e-12);
    assert(std::abs(sy - 0.02) < 1e-12);
    assert(std::abs(p.get_combined_std() - std::hypot(0.01, 0.02)) < 1e-12);

    p.set_mode(ParameterMode::SHIFTABLE);
    p.set_shift(0.05);
    assert(std::abs(p.get_val() - (0.1 + 0.05)) < 1e-12);

    p.set_std(0.2, 0.3);
    auto [st2, sy2] = p.get_std();
    assert(std::abs(st2 - 0.2) < 1e-12);
    assert(std::abs(sy2 - 0.3) < 1e-12);
}

static void test_shift_exception_on_fixed() {
    std::cout << "\n-- [Unit] shift exception on FIXED --\n";
    Parameter p(ParamId{ParameterType::SM, "GAUGE", 4}, 1.0, 0.0, 0.0);
    bool threw = false;
    try {
        p.set_shift(0.1);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw);
}

static void test_notify_observers_add_remove() {
    std::cout << "\n-- [Unit] notify observers add/remove --\n";

    auto src = std::make_shared<Parameter>(ParamId{ParameterType::SM, "SRC", 1}, 0.0, 0.0, 0.0);
    auto obs = std::make_shared<CountingParam>();

    src->addObserver(obs);
    assert(obs->updates == 0);

    src->set_expected(1.23);
    assert(obs->updates == 1);

    src->removeObserver(obs);
    src->set_expected(2.34);
    assert(obs->updates == 1);
}

static void test_assignment_and_operators() {
    std::cout << "\n-- [Unit] operator=, +=, *= and shared_ptr overloads --\n";

    Parameter a(ParamId{ParameterType::SM, "A", 10}, 1.0, 0.1, 0.2);
    a.set_mode(ParameterMode::FIXED);
    Parameter b(ParamId{ParameterType::SM, "A", 10}, 2.0, 0.3, 0.4);
    b.set_mode(ParameterMode::SHIFTABLE);
    b.set_shift(0.5);


    a = b;

    assert(std::abs(a.get_val() - (2.0 + 0.5)) < 1e-12);

    Parameter c(ParamId{ParameterType::SM, "A", 10}, 1.0, 3.0, 4.0);
    a += c;

    assert(std::abs(a.get_val() - 3.5) < 1e-12);

    a *= -2.0;

    assert(std::abs(a.get_val() - (-7.0)) < 1e-12);

    auto sp1 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "X", 1}, 1.0, 0.1, 0.2);
    auto sp2 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "X", 1}, 2.0, 0.3, 0.4);
    auto sp3 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "Y", 2}, 5.0, 0.0, 0.0);

    sp1 += sp2;
    assert(std::abs(sp1->get_val() - 3.0) < 1e-12);

    sp1 += sp3;
    assert(sp1->get_id().block == sp3->get_id().block && sp1->get_id().code == sp3->get_id().code);
    assert(std::abs(sp1->get_val() - sp3->get_val()) < 1e-12);

    sp1 * 2.0;
    assert(std::abs(sp1->get_val() - (sp3->get_val() * 2.0)) < 1e-12);
}

static void test_dependent_parameter_basic() {
    std::cout << "\n-- [Unit] DependentParameter basic --\n";

    auto p1 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "SRC", 1}, 0.5, 0.01, 0.01);
    auto p2 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "SRC", 2}, 1.0, 0.01, 0.01);

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources = {
        {p1->get_id(), p1},
        {p2->get_id(), p2}
    };

    auto lambda = [](const auto& srcs, std::shared_ptr<DependentParameter> self) {
        double sum = 0.0;
        for (const auto& [_, p] : srcs.raw()) sum += p->get_val();
        self->set_expected(sum);
    };

    auto dep = std::make_shared<DependentParameter>(ParamId{ParameterType::SM, "DEP", 100}, sources, lambda);
    dep->init();

    dep->update();
    assert(std::abs(dep->get_val() - 1.5) < 1e-12);

    p1->set_expected(1.7);
    assert(std::abs(dep->get_val() - (1.7 + 1.0)) < 1e-12);
}

static void test_dependent_freeze_unfreeze_clear_above() {
    std::cout << "\n-- [Unit] DependentParameter freeze/unfreeze/clear_above --\n";

    auto p = std::make_shared<Parameter>(ParamId{ParameterType::SM, "SRC", 1}, 1.0, 0.0, 0.0);

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> srcs = {
        {p->get_id(), p}
    };

    auto lambda = [](const auto& srcs, std::shared_ptr<DependentParameter> self) {
        self->set_expected(srcs.raw().begin()->second->get_val() * 10.0);
    };

    auto dep = std::make_shared<DependentParameter>(ParamId{ParameterType::SM, "DEP", 200}, srcs, lambda);
    dep->init();

    dep->update();
    assert(std::abs(dep->get_val() - 10.0) < 1e-12);

    dep->freeze();
    p->set_expected(2.0);
    assert(std::abs(dep->get_val() - 10.0) < 1e-12);


    dep->unfreeze();
    assert(std::abs(dep->get_val() - 20.0) < 1e-12);


    dep->clear_above();
    p->set_expected(3.0);

    assert(std::abs(dep->get_val() - 20.0) < 1e-12);
}

int main() {
    std::cout << "== Running Parameter unit tests ==\n";
    test_construct_and_std_and_mode();
    test_shift_exception_on_fixed();
    test_notify_observers_add_remove();
    test_assignment_and_operators();
    test_dependent_parameter_basic();
    test_dependent_freeze_unfreeze_clear_above();
    std::cout << "\n Parameter unit suite passed.\n";
    return 0;
}
