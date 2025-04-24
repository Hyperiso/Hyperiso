#include "DBMemento.h"
#include <iostream>
#include <cassert>

void print_block(std::shared_ptr<BlockAccessor> ba, const std::string& label) {
    std::cout << "\n[" << label << "] Block GAUGE:\n";
    for (auto& [id, val] : ba->getAllValues("GAUGE")) {
        std::cout << "  " << id << " = " << val << "\n";
    }
}

void test_dbmemento_integration_flow() {
    std::cout << "\n-- test_dbmemento_integration_flow --" << std::endl;

    DBMemento memento;

    auto ba = std::make_shared<BlockAccessor>();
    auto b = std::make_shared<Block>();
    b->blockname = "GAUGE";

    b->store(LhaID(1), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(1)), 0.1, 0, 0));
    b->store(LhaID(2), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(2)), 0.2, 0, 0));
    (*ba)["GAUGE"] = b;

    memento.takeSnapshot(ba); // Snapshot 0
    print_block(ba, "Snapshot 0");

    ba->setValue("GAUGE", LhaID(1), 1.1);
    ba->setValue("GAUGE", LhaID(2), 1.2);
    memento.takeSnapshot(ba); // Snapshot 1
    print_block(ba, "Snapshot 1");

    ba->setValue("GAUGE", LhaID(1), 2.1);
    ba->setValue("GAUGE", LhaID(2), 2.2);
    memento.takeSnapshot(ba); // Snapshot 2
    print_block(ba, "Snapshot 2");

    memento.restore(); // Restore snapshot 2 -> 1
    print_block(ba, "Restored to Snapshot 1");

    assert(std::abs(ba->getValue("GAUGE", LhaID(1)) - 1.1) < 1e-6);
    assert(std::abs(ba->getValue("GAUGE", LhaID(2)) - 1.2) < 1e-6);

    memento.restore(); // Restore snapshot 1 -> 0
    print_block(ba, "Restored to Snapshot 0");

    assert(std::abs(ba->getValue("GAUGE", LhaID(1)) - 0.1) < 1e-6);
    assert(std::abs(ba->getValue("GAUGE", LhaID(2)) - 0.2) < 1e-6);
}

int main() {
    std::cout << "== Running INTEGRATION tests for DBMemento ==\n";

    test_dbmemento_integration_flow();

    std::cout << "\n✅ All DBMemento integration tests passed!\n" << std::endl;
    return 0;
}
