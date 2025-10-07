#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "DBMemento.h"
#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "General.h"

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

int main() {
    std::cout << "== Running DBMemento integration tests ==\n";

    auto S1 = make_ba({
        {"MASS",  {{11, 100.0}, {12, 200.0}}},
        {"GAUGE", {{ 1,   0.10}}}
    });

    auto S2 = make_ba({
        {"MASS",  {{11, 101.0}, {12, 199.5}, {13, 7.7}}},
        {"GAUGE", {{ 1,   0.20}}}
    });

    auto S3 = make_ba({
        {"MASS",  {{11,  999.0}}},
        {"ALPHA", {{ 0,    0.5}}}
    });

    DBMemento mem;
    mem.takeSnapshot(S1);
    mem.takeSnapshot(S2);
    mem.takeSnapshot(S3);

    mem.restore(1);
    assert(!S3->contains("ALPHA"));
    assert(std::abs(S3->getValue("MASS", LhaID(11)) - 101.0) < 1e-12);

    assert(!S3->has_param("MASS", LhaID(13)));

    DBMemento mem2;
    mem2.takeSnapshot(S1);
    mem2.takeSnapshot(S2);
    mem2.takeSnapshot(S3);
    mem2.restore(2);

    assert(!S3->contains("ALPHA"));
    assert(std::abs(S3->getValue("MASS", LhaID(11)) - 100.0) < 1e-12);

    S3->setValue("MASS", LhaID(12), 1234.0);
    assert(std::abs(S3->getValue("MASS", LhaID(12)) - 1234.0) < 1e-12);

    S3->setValue("MASS", LhaID(14), 42.0);
    assert(S3->has_param("MASS", LhaID(14)));
    assert(std::abs(S3->getValue("MASS", LhaID(14)) - 42.0) < 1e-12);

    mem2.print_snapshot_content(0);

    std::cout << "\n✅ DBMemento integration suite passed.\n";
    return 0;
}
