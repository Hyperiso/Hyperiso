#include "Parameter.h"
#include "Block.h"
#include <iostream>
#include <cassert>
#include <cmath>

void test_Parameter_basic() {
    std::cout << "\n-- test_Parameter_basic --" << std::endl;

    ParamId pid{ParameterType::SM, "GAUGE", 3};
    Parameter p(pid, 0.1, 0.01, 0.02);

    assert(p.get_id().code == 3);
    assert(std::abs(p.get_val() - 0.1) < 1e-6);
    assert(std::abs(p.get_std() - std::hypot(0.01, 0.02)) < 1e-6);

    p.set_expected(0.2);
    assert(std::abs(p.get_val() - 0.2) < 1e-6);

    p.set_mode(ParameterMode::SHIFTABLE);
    p.set_shift(0.05);
    assert(std::abs(p.get_val() - 0.25) < 1e-6);
}

void test_DependentParameter_update() {
    std::cout << "\n-- test_DependentParameter_update --" << std::endl;

    auto p1 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "GAUGE", 1}, 0.5, 0.01, 0.01);
    auto p2 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "GAUGE", 2}, 1.0, 0.01, 0.01);

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources = {
        {p1->get_id(), p1},
        {p2->get_id(), p2}
    };

    auto lambda = [](const auto& srcs, std::shared_ptr<DependentParameter> self) {
        double sum = 0;
        for (const auto& [_, p] : srcs) sum += p->get_val();
        self->set_expected(sum);
    };

    auto dep = std::make_shared<DependentParameter>(sources, lambda);
    dep->init();
    dep->set_mode(ParameterMode::FIXED); // Just for consistent state

    p1->set_expected(0.7);  // Triggers update via notify
    assert(std::abs(dep->get_val() - 1.7) < 1e-6);
}

void test_Block_operations() {
    std::cout << "\n-- test_Block_operations --" << std::endl;

    LhaID id(5);
    auto param = std::make_shared<Parameter>();
    param->set_expected(2.5);

    Block block;
    block.blockname = "GAUGE";
    block.store(id, param);

    assert(block.contains(id));
    assert(std::abs(block.retrieve(id)->get_expected() - 2.5) < 1e-6);

    auto new_param = std::make_shared<Parameter>();
    new_param->set_expected(3.5);
    block.assign(id, new_param);
    assert(std::abs(block.retrieve(id)->get_expected() - 3.5) < 1e-6);

    block.assign(id, 4.0);
    assert(std::abs(block.retrieve(id)->get_expected() - 4.0) < 1e-6);
}

void test_DependentBlock_observer() {
    std::cout << "\n-- test_DependentBlock_observer --" << std::endl;

    LhaID id(9);
    auto source_param = std::make_shared<Parameter>();
    source_param->set_expected(10.0);

    auto source_block = std::make_shared<Block>();
    source_block->blockname = "SRC";
    source_block->store(id, source_param);

    bool triggered = false;
    auto recalc = [&](const auto& blocks, std::shared_ptr<DependentBlock> self) {
        triggered = true;
        double val = blocks.at("SRC")->retrieve(id)->get_expected();
        self->store(id, std::make_shared<Parameter>());
        self->assign(id, val + 1);
    };

    std::unordered_map<std::string, std::shared_ptr<Block>> srcs = {
        {"SRC", source_block}
    };

    auto dep = std::make_shared<DependentBlock>(srcs, recalc);
    dep->blockname = "DEP";
    dep->init();

    source_block->assign(id, 15.0);
    assert(triggered);
    assert(dep->contains(id));
    assert(std::abs(dep->retrieve(id)->get_expected() - 16.0) < 1e-6);
}

void test_WilsonBlock_access() {
    std::cout << "\n-- test_WilsonBlock_access --" << std::endl;

    WilsonBlock wb;
    LhaID coefID(30);  // 3 * 10 + 0
    wb.setValue(coefID, 0.25);
    assert(std::abs(wb.getValue(coefID) - 0.25) < 1e-6);

    wb.setValue(LhaID(-1), 1000);  // Scale
    wb.setValue(LhaID(-2), 42);    // Type
    assert(std::abs(wb.getValue(LhaID(-1)) - 1000) < 1e-6);
    assert(std::abs(wb.getValue(LhaID(-2)) - 42) < 1e-6);
}

int main() {
    std::cout << "== Launching full unit test for Parameter & Block system ==" << std::endl;

    test_Parameter_basic();
    test_DependentParameter_update();
    test_Block_operations();
    test_DependentBlock_observer();
    test_WilsonBlock_access();

    std::cout << "\n✅ All tests passed with success!\n" << std::endl;
    return 0;
}
