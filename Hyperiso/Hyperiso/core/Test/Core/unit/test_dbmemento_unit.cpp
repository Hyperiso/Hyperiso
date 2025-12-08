#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "DBMemento.h"
#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "Include.h"

static std::shared_ptr<Block> make_block(const BlockName& name,
                                         const std::vector<std::pair<long,double>>& vals) {
    auto b = std::make_shared<Block>();
    b->blockname = name;
    for (auto& [id, v] : vals) {
        auto p = std::make_shared<Parameter>(ParamId{ParameterType::SM, name, id}, v, 0.0, 0.0);
        b->store(LhaID(id), p);
    }
    return b;
}

static std::shared_ptr<BlockAccessor> make_ba(
    const std::vector<std::pair<BlockName, std::vector<std::pair<long,double>>>>& spec)
{
    auto ba = std::make_shared<BlockAccessor>();
    for (auto& [bn, vec] : spec) {
        auto blk = make_block(bn, vec);
        ba->emplace(bn, blk);
    }
    return ba;
}

static void test_stack_and_restore_1step() {
    std::cout << "\n-- unit: stack + restore(1) --\n";

    auto ba1 = make_ba({{"MASS",  {{11,1.0}}},
                        {"GAUGE", {{ 1,0.5}}}});
    auto ba2 = make_ba({{"MASS",  {{11,2.0}, {12, 7.}}}, 
                        {"GAUGE", {{ 1,0.7}}}}); 
    auto ba3 = make_ba({{"MASS",  {{11,9.9}}},   
                        {"ALPHA", {{ 0,0.1}}}});

    DBMemento mem;
    assert(mem.stack_size() == 0);
    mem.takeSnapshot(ba1);
    mem.takeSnapshot(ba2);
    mem.takeSnapshot(ba3);
    assert(mem.stack_size() == 3);

    mem.restore(1);
    assert(mem.stack_size() == 2); 

    assert(!ba3->contains("ALPHA"));
    assert( std::abs(ba3->getValue("MASS", LhaID(11)) - 2.0) < 1e-12 );

    bool has_mass12 = ba3->has_param("MASS", LhaID(12));
    assert(!has_mass12);
}

static void test_restore_multi_steps() {
    std::cout << "\n-- unit: restore(2) multi-steps --\n";
    auto baA = make_ba({{"GAUGE", {{1, 0.1}}},
                        {"MASS",  {{11, 100.0}}}});
    auto baB = make_ba({{"GAUGE", {{1, 0.2}}},
                        {"MASS",  {{11, 110.0}}}});
    auto baC = make_ba({{"GAUGE", {{1, 0.9}}},
                        {"MASS",  {{11, 999.0}}},
                        {"EXTRA", {{5,  1.23}}}});

    DBMemento mem;
    mem.takeSnapshot(baA);
    mem.takeSnapshot(baB);
    mem.takeSnapshot(baC);

    mem.restore(2);

    assert(!baC->contains("EXTRA"));

    assert(std::abs(baC->getValue("GAUGE", LhaID(1)) - 0.1) < 1e-12);
    assert(std::abs(baC->getValue("MASS",  LhaID(11)) - 100.0) < 1e-12);
}

static void test_print_snapshot_content_no_crash() {
    std::cout << "\n-- unit: print_snapshot_content no-crash --\n";
    auto ba1 = make_ba({{"MASS", {{11, 1.0}}}});
    auto ba2 = make_ba({{"MASS", {{11, 2.0}}}});

    DBMemento mem;
    mem.takeSnapshot(ba1);
    mem.takeSnapshot(ba2);

    mem.print_snapshot_content(0);
    mem.print_snapshot_content(1);
}

int main() {
    std::cout << "== Running DBMemento unit tests ==\n";
    test_stack_and_restore_1step();
    test_restore_multi_steps();
    test_print_snapshot_content_no_crash();
    std::cout << "\n DBMemento unit suite passed.\n";
    return 0;
}
