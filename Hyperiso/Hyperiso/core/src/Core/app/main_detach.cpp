#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "BlockAccessor.h"
#include "Block.h"
#include "DependentParameter.h"
#include "SourcesView.h"

namespace {

using PT = ParameterType;

bool almost_equal(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) < eps;
}

void check(BlockAccessor& ba, const std::string& block, const LhaID& id, double expected) {
    const double got = ba.getValue(block, id);
    std::cout << block << "[" << id << "] = " << got << " (expected " << expected << ")\n";
    assert(almost_equal(got, expected));
}

void banner(const std::string& title) {
    std::cout << "\n==============================\n";
    std::cout << title << "\n";
    std::cout << "==============================\n";
}

void add_leaf_block(BlockAccessor& ba, const BlockName& name) {
    auto blk = std::make_shared<Block>();
    blk->blockname = name;
    ba.emplace(name, blk);
}

void set_typed_leaf(BlockAccessor& ba, const BlockName& block, const LhaID& id, double value, PT type) {
    ba.setValue(block, id, value);
    ba.getParameter(block, id)->set_owner(type);
}

void add_dependent_block_local(
    BlockAccessor& ba,
    const BlockName& name,
    std::unordered_map<std::string, std::shared_ptr<Block>> sources,
    DepUpdateFunc recalculateFunc)
{
    auto dep = std::make_shared<DependentBlock>(sources, std::move(recalculateFunc));
    dep->blockname = name;
    ba.emplace(name, dep);   // important avant init()
    dep->init();
    dep->update();
}

void add_dependent_parameter_local(
    BlockAccessor& ba,
    ParamId pid,
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources,
    DepParamUpdateFunc recalculateFunc)
{
    if (!ba.contains(pid.block)) {
        add_leaf_block(ba, pid.block);
    }

    auto blk = ba.at(pid.block);

    if (blk->contains(pid.code)) {
        auto existing = blk->retrieve(pid.code);
        if (auto dep = std::dynamic_pointer_cast<DependentParameter>(existing)) {
            dep->rebind(std::move(sources), std::move(recalculateFunc));
            dep->init();
            return;
        }
    }

    auto dep = std::make_shared<DependentParameter>(pid, std::move(sources), std::move(recalculateFunc));
    blk->store(pid.code, dep);
    dep->init();
    dep->update();
}

} // namespace

int main() {
    BlockAccessor ba;

    const LhaID A{1};
    const LhaID C{1};
    const LhaID P{1};
    const LhaID Z{1};

    // -------------------------------------------------------------------------
    // 1) BASE : block leaf
    // -------------------------------------------------------------------------
    add_leaf_block(ba, "BASE");
    set_typed_leaf(ba, "BASE", A, 10.0, PT::SM);

    // -------------------------------------------------------------------------
    // 2) DB1 : DependentBlock <- BASE
    //
    // c = 2 * BASE[a]
    // -------------------------------------------------------------------------
    add_dependent_block_local(
        ba,
        "DB1",
        {
            {"BASE", ba.at("BASE")}
        },
        [](const BlockSrc& src, std::shared_ptr<DependentBlock> self) {
            const double a = src.get_val("BASE", 1);
            const double c = 2.0 * a;

            auto out = std::make_shared<Parameter>(
                ParamId(PT::SM, "DB1", LhaID{1}),
                c, 0.0, 0.0
            );
            self->store_or_assign(LhaID{1}, out);
        }
    );

    // force existence de DB1[1]
    (void)ba.getValue("DB1", C);

    // -------------------------------------------------------------------------
    // 3) PARAMS : block leaf contenant un DependentParameter
    //
    // p = DB1[c] + 5
    // -------------------------------------------------------------------------
    {
        ParamId src_id(PT::SM, "DB1", C);
        ParamId dst_id(PT::SM, "PARAMS", P);

        add_dependent_parameter_local(
            ba,
            dst_id,
            {
                {src_id, ba.getParameter("DB1", C)}
            },
            [src_id](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
                const double c = src.get_val(src_id);
                self->set_expected_silent(c + 5.0);
            }
        );
    }

    // -------------------------------------------------------------------------
    // 4) DB2 : DependentBlock <- DB1 + PARAMS
    //
    // z = DB1[c] + PARAMS[p]
    // -------------------------------------------------------------------------
    add_dependent_block_local(
        ba,
        "DB2",
        {
            {"DB1", ba.at("DB1")},
            {"PARAMS", ba.at("PARAMS")}
        },
        [](const BlockSrc& src, std::shared_ptr<DependentBlock> self) {
            const double c = src.get_val("DB1", 1);
            const double p = src.get_val("PARAMS", 1);
            const double z = c + p;

            auto out = std::make_shared<Parameter>(
                ParamId(PT::SM, "DB2", LhaID{1}),
                z, 0.0, 0.0
            );
            self->store_or_assign(LhaID{1}, out);
        }
    );

    // -------------------------------------------------------------------------
    // Etat initial
    // BASE[a] = 10
    // DB1[c] = 20
    // PARAMS[p] = 25
    // DB2[z] = 45
    // -------------------------------------------------------------------------
    banner("Etat initial");
    check(ba, "BASE",   A, 10.0);
    check(ba, "DB1",    C, 20.0);
    check(ba, "PARAMS", P, 25.0);
    check(ba, "DB2",    Z, 45.0);

    // -------------------------------------------------------------------------
    // BASE[a] = 20
    // DB1[c] = 40
    // PARAMS[p] = 45
    // DB2[z] = 85
    // -------------------------------------------------------------------------
    banner("Apres BASE[a] = 20");
    ba.setValue("BASE", A, 20.0);
    check(ba, "BASE",   A, 20.0);
    check(ba, "DB1",    C, 40.0);
    check(ba, "PARAMS", P, 45.0);
    check(ba, "DB2",    Z, 85.0);

    // -------------------------------------------------------------------------
    // detach PARAMS[p], puis BASE[a] = 7
    // DB1[c] = 14
    // PARAMS[p] reste 45
    // DB2[z] = 59
    // -------------------------------------------------------------------------
    banner("Detach PARAMS[p], puis BASE[a] = 7");
    ba.detach_parameter("PARAMS", P);
    ba.setValue("BASE", A, 7.0);
    check(ba, "BASE",   A, 7.0);
    check(ba, "DB1",    C, 14.0);
    check(ba, "PARAMS", P, 45.0);
    check(ba, "DB2",    Z, 59.0);

    // -------------------------------------------------------------------------
    // reattach PARAMS[p]
    // DB1[c] = 14
    // PARAMS[p] = 19
    // DB2[z] = 33
    // -------------------------------------------------------------------------
    banner("Reattach PARAMS[p]");
    ba.reattach_parameter("PARAMS", P);
    check(ba, "BASE",   A, 7.0);
    check(ba, "DB1",    C, 14.0);
    check(ba, "PARAMS", P, 19.0);
    check(ba, "DB2",    Z, 33.0);

    // -------------------------------------------------------------------------
    // detach DB1, puis BASE[a] = 100
    // DB1[c] reste 14
    // PARAMS[p] reste 19
    // DB2[z] reste 33
    // -------------------------------------------------------------------------
    banner("Detach DB1, puis BASE[a] = 100");
    ba.detach_block("DB1");
    ba.setValue("BASE", A, 100.0);
    check(ba, "BASE",   A, 100.0);
    check(ba, "DB1",    C, 14.0);
    check(ba, "PARAMS", P, 19.0);
    check(ba, "DB2",    Z, 33.0);

    // -------------------------------------------------------------------------
    // reattach DB1
    // DB1[c] = 200
    // PARAMS[p] = 205
    // DB2[z] = 405
    // -------------------------------------------------------------------------
    banner("Reattach DB1");
    ba.reattach_block("DB1");
    check(ba, "BASE",   A, 100.0);
    check(ba, "DB1",    C, 200.0);
    check(ba, "PARAMS", P, 205.0);
    check(ba, "DB2",    Z, 405.0);

    std::cout << "\nTous les checks sont passes.\n";
    return 0;
}