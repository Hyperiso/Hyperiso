#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Block.h"
#include "Parameter.h"
#include "Include.h"
#include "SourcesView.h"

struct CountingParameter : Parameter {
    int updates = 0;
    int freezes = 0;
    int unfreezes = 0;
    int clear_above_calls = 0;
    int clear_below_calls = 0;

    CountingParameter(const ParamId& id, double v)
        : Parameter(id, v, 0.0, 0.0) {}

    void update() override { ++updates; }
    void freeze() override { ++freezes; }
    void unfreeze() override { ++unfreezes; }
    void clear_above() override { ++clear_above_calls; }
    void clear_below() override { ++clear_below_calls; }
};

struct CountingBlock : Block {
    int updates = 0;
    void update() override { ++updates; }
};

static std::shared_ptr<Parameter> mkp(const std::string& blk, int code, double v) {
    return std::make_shared<Parameter>(ParamId{ParameterType::SM, blk, code}, v, 0.0, 0.0);
}

static void test_basic_store_retrieve_assign_remove_contains() {
    std::cout << "\n-- [Unit] Block: store/retrieve/assign/remove/contains --\n";

    auto b = std::make_shared<Block>();
    b->blockname = "GAUGE";
    b->bind_self(b);

    LhaID id(11);
    auto p = std::make_shared<Parameter>(ParamId{ParameterType::SM, "GAUGE", 11}, 1.0, 0.0, 0.0);

    b->store(id, p);
    assert(b->contains(id));
    assert(std::abs(b->retrieve(id)->get_val() - 1.0) < 1e-12);

    auto p2 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "GAUGE", 11}, 2.5, 0.0, 0.0);
    b->assign(id, p2);
    assert(std::abs(b->retrieve(id)->get_val() - 2.5) < 1e-12);

    b->assign(id, 3.25);
    assert(std::abs(b->retrieve(id)->get_val() - 3.25) < 1e-12);

    LhaID id2(12);
    auto pNew = std::make_shared<Parameter>(ParamId{ParameterType::SM, "GAUGE", 12}, 7.0, 0.0, 0.0);
    b->store_or_assign(id2, pNew);
    assert(b->contains(id2));
    assert(std::abs(b->retrieve(id2)->get_val() - 7.0) < 1e-12);

    b->remove(id2);
    assert(!b->contains(id2));
}

static void test_getAllIDs_and_getItems() {
    std::cout << "\n-- [Unit] Block: getAllIDs/getItems --\n";

    auto b = std::make_shared<Block>();
    b->blockname = "TEST";
    b->bind_self(b);

    auto p1 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "TEST", 1}, 1.0, 0.0, 0.0);
    auto p2 = std::make_shared<Parameter>(ParamId{ParameterType::SM, "TEST", 2}, 2.0, 0.0, 0.0);

    b->store(LhaID(1), p1);
    b->store(LhaID(2), p2);

    auto ids = b->getAllIDs();
    assert(ids.count(LhaID(1)) == 1);
    assert(ids.count(LhaID(2)) == 1);
    assert(ids.size() == 2);

    const auto& items = b->getItems();
    assert(items.size() == 2);
}

static void test_update_freeze_unfreeze_delegate_to_params() {
    std::cout << "\n-- [Unit] Block: update/freeze/unfreeze delegates --\n";

    auto b = std::make_shared<Block>();
    b->blockname = "DELEG";
    b->bind_self(b);

    auto c1 = std::make_shared<CountingParameter>(ParamId{ParameterType::SM, "DELEG", 1}, 0.0);
    auto c2 = std::make_shared<CountingParameter>(ParamId{ParameterType::SM, "DELEG", 2}, 0.0);
    b->store(LhaID(1), c1);
    b->store(LhaID(2), c2);

    b->update();
    assert(c1->updates == 1 && c2->updates == 1);

    b->freeze();
    assert(c1->freezes == 1 && c2->freezes == 1);

    b->unfreeze();
    assert(c1->unfreezes == 1 && c2->unfreezes == 1);
}

static void test_observers_notify_between_blocks() {
    std::cout << "\n-- [Unit] Block: notifyObservers to other Block --\n";

    auto src = std::make_shared<Block>();
    src->blockname = "SRC";
    src->bind_self(src);

    auto dst = std::make_shared<CountingBlock>();
    dst->blockname = "DST";
    dst->bind_self(dst);

    src->addObserver(dst);
    assert(dst->updates == 0);

    LhaID k(42);
    src->store(k, std::make_shared<Parameter>(
        ParamId{ParameterType::SM, "SRC", 42}, 0.0, 0.0, 0.0
    ));

    src->assign(k, 3.14);
    assert(dst->updates == 1);

    src->removeObserver(dst);
    src->assign(k, 2.71);
    assert(dst->updates == 1);
}

static void test_copy_ctor_and_copy_method_and_scale() {
    std::cout << "\n-- [Unit] Block: copy/ctor + scale --\n";

    auto src = std::make_shared<Block>();
    src->blockname = "SRC";
    src->bind_self(src);

    src->store(LhaID(5), std::make_shared<Parameter>(ParamId{ParameterType::SM, "SRC", 5}, 9.0, 0.0, 0.0));
    assert(!src->has_scale());
    src->set_scale(100.0);
    assert(src->has_scale());
    assert(std::abs(src->get_scale() - 100.0) < 1e-12);

    auto dst = std::make_shared<Block>();
    dst->bind_self(dst);
    dst->copy(src);

    assert(dst->contains(LhaID(5)));
    assert(std::abs(dst->retrieve(LhaID(5))->get_val() - 9.0) < 1e-12);
    assert(dst->has_scale());
    assert(std::abs(dst->get_scale() - 100.0) < 1e-12);

    Block cpy(src);
    assert(cpy.contains(LhaID(5)));
    assert(std::abs(cpy.retrieve(LhaID(5))->get_val() - 9.0) < 1e-12);
    assert(cpy.has_scale());
    assert(std::abs(cpy.get_scale() - 100.0) < 1e-12);
}

static void test_clear_above_below() {
    std::cout << "\n-- [Unit] Block: clear_above/clear_below --\n";

    auto b = std::make_shared<Block>();
    b->blockname = "CLR";
    b->bind_self(b);

    auto c1 = std::make_shared<CountingParameter>(ParamId{ParameterType::SM, "CLR", 1}, 0.0);
    auto c2 = std::make_shared<CountingParameter>(ParamId{ParameterType::SM, "CLR", 2}, 0.0);
    b->store(LhaID(1), c1);
    b->store(LhaID(2), c2);

    b->clear_above();
    assert(c1->clear_above_calls == 1 && c2->clear_above_calls == 1);

    b->clear_below();
    assert(c1->clear_below_calls == 1 && c2->clear_below_calls == 1);
}

static void test_dependent_block_lifecycle_lazy_dirty_freeze_clear_above() {
    std::cout << "\n-- [Unit] DependentBlock: lazy/dirty/freeze/unfreeze/clear_above --\n";

    auto src = std::make_shared<Block>();
    src->blockname = "SRC";
    src->bind_self(src);

    LhaID key(7);
    src->store(key, std::make_shared<Parameter>(ParamId{ParameterType::SM, "SRC", 7}, 1.0, 0.0, 0.0));

    int runs = 0;
    auto recalc = [&](const auto& blocks, std::shared_ptr<DependentBlock> self) {
        ++runs;
        double v = blocks.get_val("SRC", key);
        if (!self->contains(key)) self->store(key, std::make_shared<Parameter>(ParamId{ParameterType::SM, "DEP", 7}, 0.0, 0.0, 0.0));
        self->assign(key, v + 10.0);
    };

    std::unordered_map<std::string, std::shared_ptr<Block>> sources = { {"SRC", src} };
    auto dep = std::make_shared<DependentBlock>(sources, recalc);
    dep->blockname = "DEP";
    dep->bind_self(dep);
    dep->init();

    assert(dep->dependsOn("SRC"));
    assert(!dep->dependsOn("NOPE"));

    assert(runs == 0);

    dep->freeze();
    src->assign(key, 5.0);
    assert(runs == 0);

    dep->unfreeze();
    assert(runs == 0);

    double v = dep->retrieve(key)->get_val();
    assert(runs == 1);
    assert(std::abs(v - 15.0) < 1e-12);

    dep->clear_above();
    src->assign(key, 9.0);

    double v2 = dep->retrieve(key)->get_val();
    assert(runs == 1);
    assert(std::abs(v2 - 15.0) < 1e-12);
}

static void test_dependent_block_detach_reattach() {
    std::cout << "\n-- [Unit] DependentBlock detach/reattach --\n";

    auto src = std::make_shared<Block>();
    src->blockname = "SRCD";
    src->bind_self(src);

    LhaID key(8);
    src->store(key, std::make_shared<Parameter>(
        ParamId{ParameterType::SM, "SRCD", 8}, 2.0, 0.0, 0.0
    ));

    int runs = 0;
    auto recalc = [&](const auto& blocks, std::shared_ptr<DependentBlock> self) {
        ++runs;
        double v = blocks.get_val("SRCD", key);
        if (!self->contains(key)) {
            self->store(key, std::make_shared<Parameter>(
                ParamId{ParameterType::SM, "DETD", 8}, 0.0, 0.0, 0.0
            ));
        }
        self->assign(key, v + 100.0);
    };

    auto dep = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{{"SRCD", src}},
        recalc
    );
    dep->blockname = "DETD";
    dep->bind_self(dep);
    dep->init();

    assert(std::abs(dep->retrieve(key)->get_val() - 102.0) < 1e-12);
    assert(runs == 1);

    dep->detach();

    src->assign(key, 5.0);

    assert(std::abs(dep->retrieve(key)->get_val() - 102.0) < 1e-12);
    assert(runs == 1);

    dep->reattach();

    assert(std::abs(dep->retrieve(key)->get_val() - 105.0) < 1e-12);
    assert(runs == 2);

    src->assign(key, 9.0);
    assert(std::abs(dep->retrieve(key)->get_val() - 109.0) < 1e-12);
    assert(runs == 3);
}

int main() {
    std::cout << "== Running Block unit tests ==\n";
    test_basic_store_retrieve_assign_remove_contains();
    test_getAllIDs_and_getItems();
    test_update_freeze_unfreeze_delegate_to_params();
    test_observers_notify_between_blocks();
    test_copy_ctor_and_copy_method_and_scale();
    test_clear_above_below();
    test_dependent_block_lifecycle_lazy_dirty_freeze_clear_above();
    test_dependent_block_detach_reattach();
    std::cout << "\n Block unit suite passed.\n";
    return 0;
}
