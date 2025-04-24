#include "BlockAccessor.h"
#include "Parameter.h"
#include "General.h"
#include <iostream>
#include <cassert>
#include <cmath>

void test_full_parameter_flow() {
    std::cout << "\n-- test_full_parameter_flow --" << std::endl;

    BlockAccessor accessor;

    std::string block_name = "GAUGE";
    LhaID id(3);  // Could represent g3
    ParamId pid(ParameterType::SM, block_name, id);

    auto param = std::make_shared<Parameter>(pid, 0.118, 0.001, 0.002);
    auto block = std::make_shared<Block>();
    block->blockname = block_name;

    // Store and attach block to accessor
    block->store(id, param);
    accessor[block_name] = block;

    // Check stored value
    assert(accessor.has_param(block_name, id));
    auto retrieved = accessor.getParameter(block_name, id);
    assert(std::abs(retrieved->get_val() - 0.118) < 1e-6);

    // Assign new value
    accessor.setValue(block_name, id, 0.120);
    assert(std::abs(accessor.getValue(block_name, id) - 0.120) < 1e-6);

    // Remove
    accessor.remove_item(block_name, id);
    assert(!accessor.has_param(block_name, id));
}

void test_multikey_lhaid() {
    std::cout << "\n-- test_multikey_lhaid --" << std::endl;

    BlockAccessor accessor;
    std::string block_name = "MATRIX";
    LhaID id(1, 2);  // simulate a matrix entry

    ParamId pid(ParameterType::SM, block_name, id);
    auto block = std::make_shared<Block>();
    block->blockname = block_name;

    auto param = std::make_shared<Parameter>(pid, 42.0, 0.1, 0.2);
    block->store(id, param);
    accessor[{block_name}] = block;

    assert(accessor.has_param(block_name, id));
    assert(std::abs(accessor.getValue(block_name, id) - 42.0) < 1e-6);
}

void test_getAllValues_and_subaccessor() {
    std::cout << "\n-- test_getAllValues_and_subaccessor --" << std::endl;

    BlockAccessor accessor;

    // Block 1
    auto b1 = std::make_shared<Block>();
    b1->blockname = "BLOCK1";
    b1->store(LhaID(1), std::make_shared<Parameter>(ParamId("BLOCK1", LhaID(1)), 10.0, 0, 0));
    b1->store(LhaID(2), std::make_shared<Parameter>(ParamId("BLOCK1", LhaID(2)), 20.0, 0, 0));
    accessor["BLOCK1"] = b1;

    // Block 2
    auto b2 = std::make_shared<Block>();
    b2->blockname = "BLOCK2";
    b2->store(LhaID(1), std::make_shared<Parameter>(ParamId("BLOCK2", LhaID(1)), 30.0, 0, 0));
    accessor["BLOCK2"] = b2;

    auto values = accessor.getAllValues("BLOCK1");
    assert(values.size() == 2);
    assert(std::abs(values[LhaID(1)] - 10.0) < 1e-6);
    assert(std::abs(values[LhaID(2)] - 20.0) < 1e-6);

    // Subset
    auto subset = accessor[{ "BLOCK2" }];
    assert(subset->size() == 1);
    assert(subset->contains("BLOCK2"));
    assert(std::abs(subset->getValue("BLOCK2", LhaID(1)) - 30.0) < 1e-6);
}

void test_merge_behaviors() {
    std::cout << "\n-- test_merge_behaviors --" << std::endl;

    auto lhs = std::make_shared<BlockAccessor>();
    auto rhs = std::make_shared<BlockAccessor>();

    auto block1 = std::make_shared<Block>();
    block1->blockname = "GAUGE";
    block1->store(LhaID(1), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(1)), 0.1, 0, 0));
    (*lhs)["GAUGE"] = block1;

    auto block2 = std::make_shared<Block>();
    block2->blockname = "GAUGE";
    block2->store(LhaID(1), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(1)), 0.2, 0, 0));
    (*rhs)["GAUGE"] = block2;

    auto merged = *lhs >> rhs;
    assert(std::abs(merged->getValue("GAUGE", LhaID(1)) - 0.1) < 1e-6);

    try {
        auto conflict = lhs + rhs;
        assert(false && "Expected error on + merge");
    } catch (const std::exception& e) {
        std::cout << "Caught expected error on merge + : " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "== Running INTEGRATION tests for BlockAccessor ==\n";

    test_full_parameter_flow();
    test_multikey_lhaid();
    test_getAllValues_and_subaccessor();
    test_merge_behaviors();

    std::cout << "\n✅ All BlockAccessor integration tests passed!\n" << std::endl;
    return 0;
}
