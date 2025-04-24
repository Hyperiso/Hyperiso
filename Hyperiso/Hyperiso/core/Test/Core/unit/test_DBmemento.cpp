#include "DBMemento.h"
#include <iostream>
#include <cassert>
#include <memory>

// Utilitaire simplifié pour créer un BlockAccessor
std::shared_ptr<BlockAccessor> make_block_accessor(double val1, double val2) {
    auto ba = std::make_shared<BlockAccessor>();

    auto b = std::make_shared<Block>();
    b->blockname = "GAUGE";
    b->store(LhaID(1), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(1)), val1, 0, 0));
    b->store(LhaID(2), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(2)), val2, 0, 0));

    (*ba)["GAUGE"] = b;
    return ba;
}

void test_takeSnapshot_and_stack_size() {
    std::cout << "\n-- test_takeSnapshot_and_stack_size --" << std::endl;

    DBMemento mem;

    auto ba1 = make_block_accessor(1.0, 2.0);
    auto ba2 = make_block_accessor(3.0, 4.0);

    mem.takeSnapshot(ba1);
    mem.takeSnapshot(ba2);

    assert(mem.stack_size() == 2);
}

void test_restore_overwrites_values() {
    std::cout << "\n-- test_restore_overwrites_values --" << std::endl;

    DBMemento mem;

    auto ba1 = make_block_accessor(1.0, 2.0);
    auto ba2 = make_block_accessor(10.0, 20.0); // valeurs différentes

    mem.takeSnapshot(ba1);
    mem.takeSnapshot(ba2); // HEAD

    mem.restore(); // restore ba2 -> ba1

    auto restored = ba2;
    assert(std::abs(restored->getValue("GAUGE", LhaID(1)) - 1.0) < 1e-6);
    assert(std::abs(restored->getValue("GAUGE", LhaID(2)) - 2.0) < 1e-6);
}

void test_restore_multiple_steps() {
    std::cout << "\n-- test_restore_multiple_steps --" << std::endl;

    DBMemento mem;
    auto ba1 = make_block_accessor(1.0, 1.0);
    auto ba2 = make_block_accessor(2.0, 2.0);
    auto ba3 = make_block_accessor(3.0, 3.0);

    mem.takeSnapshot(ba1);
    mem.takeSnapshot(ba2);
    mem.takeSnapshot(ba3); // HEAD

    mem.restore(2); // ba3 sera remplacé par ba1

    auto restored = ba3;
    assert(std::abs(restored->getValue("GAUGE", LhaID(1)) - 1.0) < 1e-6);
    assert(std::abs(restored->getValue("GAUGE", LhaID(2)) - 1.0) < 1e-6);
}

int main() {
    std::cout << "== Running UNIT tests for DBMemento ==\n";

    test_takeSnapshot_and_stack_size();
    test_restore_overwrites_values();
    test_restore_multiple_steps();

    std::cout << "\n✅ All DBMemento unit tests passed!\n" << std::endl;
    return 0;
}
