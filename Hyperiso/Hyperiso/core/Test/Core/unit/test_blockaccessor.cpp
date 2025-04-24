#include "BlockAccessor.h"
#include <iostream>
#include <cassert>
#include <memory>
#include <cmath>

// === MOCK PARAMETER (simplié, mais pas par héritage) ===

class MockParameter : public Parameter {
public:
    MockParameter(double val) {
        set_expected(val);
    }
};

// === UNIT TESTS ===

void test_add_and_getValue() {
    std::cout << "\n-- test_add_and_getValue --" << std::endl;

    BlockAccessor ba;

    auto block = std::make_shared<Block>();
    block->blockname = "GAUGE";

    LhaID id(1);
    auto param = std::make_shared<MockParameter>(0.118);
    block->store(id, param);
    ba["GAUGE"] = block;

    assert(ba.has_param("GAUGE", id));
    assert(std::abs(ba.getValue("GAUGE", id) - 0.118) < 1e-6);
}

void test_setValue_and_getParameter() {
    std::cout << "\n-- test_setValue_and_getParameter --" << std::endl;

    BlockAccessor ba;
    auto block = std::make_shared<Block>();
    block->blockname = "FERMION";

    LhaID id(6);
    ba[{"FERMION"}] = block;

    ba.setValue("FERMION", id, 173.0);
    assert(ba.has_param("FERMION", id));
    assert(std::abs(ba.getValue("FERMION", id) - 173.0) < 1e-6);

    auto p = ba.getParameter("FERMION", id);
    assert(std::abs(p->get_val() - 173.0) < 1e-6);
}

void test_getAllValues_and_get_block_names() {
    std::cout << "\n-- test_getAllValues_and_get_block_names --" << std::endl;

    BlockAccessor ba;
    auto block = std::make_shared<Block>();
    block->blockname = "SMINPUTS";

    LhaID id1(1), id2(2);
    block->store(id1, std::make_shared<MockParameter>(0.007297));
    block->store(id2, std::make_shared<MockParameter>(91.1876));
    ba["SMINPUTS"] = block;

    auto values = ba.getAllValues("SMINPUTS");
    assert(values.size() == 2);
    assert(std::abs(values[id1] - 0.007297) < 1e-6);
    assert(std::abs(values[id2] - 91.1876) < 1e-6);

    auto names = ba.get_block_names();
    assert(names.contains("SMINPUTS"));
}

void test_remove_item() {
    std::cout << "\n-- test_remove_item --" << std::endl;

    BlockAccessor ba;
    auto block = std::make_shared<Block>();
    block->blockname = "TO_REMOVE";

    LhaID id(10);
    block->store(id, std::make_shared<MockParameter>(3.14));
    ba["TO_REMOVE"] = block;

    ba.remove_item("TO_REMOVE", id);
    assert(!ba.has_param("TO_REMOVE", id));
}

void test_operator_plus_error() {
    std::cout << "\n-- test_operator_plus_error --" << std::endl;

    auto ba1 = std::make_shared<BlockAccessor>();
    auto ba2 = std::make_shared<BlockAccessor>();

    auto block1 = std::make_shared<Block>();
    block1->blockname = "DUPLICATE";
    block1->store(LhaID(1), std::make_shared<MockParameter>(1.0));
    (*ba1)["DUPLICATE"] = block1;

    auto block2 = std::make_shared<Block>();
    block2->blockname = "DUPLICATE";
    block2->store(LhaID(2), std::make_shared<MockParameter>(2.0));
    (*ba2)["DUPLICATE"] = block2;

    try {
        auto merged = ba1 + ba2;
        assert(false && "Expected error for duplicate block");
    } catch (const std::exception& e) {
        std::cout << "Caught expected error: " << e.what() << std::endl;
    }
}

void test_operator_shift_merge() {
    std::cout << "\n-- test_operator_shift_merge --" << std::endl;

    auto ba1 = std::make_shared<BlockAccessor>();
    auto ba2 = std::make_shared<BlockAccessor>();

    auto b1 = std::make_shared<Block>();
    b1->blockname = "GAUGE";
    b1->store(LhaID(1), std::make_shared<MockParameter>(0.1));
    (*ba1)["GAUGE"] = b1;

    auto b2 = std::make_shared<Block>();
    b2->blockname = "GAUGE";
    b2->store(LhaID(1), std::make_shared<MockParameter>(0.2));
    (*ba2)["GAUGE"] = b2;

    auto merged = ba1 >> ba2;
    assert(std::abs(merged->getValue("GAUGE", LhaID(1)) - 0.1) < 1e-6);
}

void test_operator_bracket_subaccessor() {
    std::cout << "\n-- test_operator_bracket_subaccessor --" << std::endl;

    BlockAccessor ba;

    auto b1 = std::make_shared<Block>();
    b1->blockname = "BLOCK1";
    ba["BLOCK1"] = b1;

    auto b2 = std::make_shared<Block>();
    b2->blockname = "BLOCK2";
    ba["BLOCK2"] = b2;

    std::unordered_set<std::string> names = {"BLOCK1"};
    auto sub = ba[names];

    assert(sub->size() == 1);
    assert(sub->contains("BLOCK1"));
}

int main() {
    std::cout << "== Running UNIT tests for BlockAccessor (with actual Block, mock param) ==\n";

    test_add_and_getValue();
    test_setValue_and_getParameter();
    test_getAllValues_and_get_block_names();
    test_remove_item();
    test_operator_plus_error();
    test_operator_shift_merge();
    test_operator_bracket_subaccessor();

    std::cout << "\n✅ All BlockAccessor unit tests passed!\n" << std::endl;
    return 0;
}
