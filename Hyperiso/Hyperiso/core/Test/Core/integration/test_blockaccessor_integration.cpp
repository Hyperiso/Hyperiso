#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "Include.h"
#include "SourcesView.h"
#include "DependentParameter.h"

static void add_leaf_block(BlockAccessor& ba, const BlockName& name) {
    auto blk = std::make_shared<Block>();
    blk->blockname = name;
    ba.emplace(name, blk);
}

static void set_typed_leaf(BlockAccessor& ba, const BlockName& block, const LhaID& id, double value) {
    ba.setValue(block, id, value);
    ba.getParameter(block, id)->set_owner(ParameterType::SM);
}

static void add_dependent_block_local(
    BlockAccessor& ba,
    const BlockName& name,
    std::unordered_map<std::string, std::shared_ptr<Block>> sources,
    DepUpdateFunc recalculateFunc)
{
    auto dep = std::make_shared<DependentBlock>(sources, std::move(recalculateFunc));
    dep->blockname = name;
    ba.emplace(name, dep);
    dep->init();
    dep->update();
}

static void add_dependent_parameter_local(
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

static std::shared_ptr<Parameter> mkp(const std::string& blk, long code, double v){
    return std::make_shared<Parameter>(ParamId{ParameterType::SM, blk, LhaID(code)}, v, 0.0, 0.0);
}

static void test_block_accessor_mixed_chain_detach_reattach() {
    std::cout << "\n-- [Integration] BlockAccessor mixed chain detach/reattach --\n";

    BlockAccessor ba;

    const LhaID A{1};
    const LhaID C{1};
    const LhaID P{1};
    const LhaID Z{1};

    add_leaf_block(ba, "BASE");
    set_typed_leaf(ba, "BASE", A, 10.0);

    // DB1[c] = 2 * BASE[a]
    add_dependent_block_local(
        ba,
        "DB1",
        {
            {"BASE", ba.at("BASE")}
        },
        [](const BlockSrc& src, std::shared_ptr<DependentBlock> self) {
            double a = src.get_val("BASE", 1);
            auto out = std::make_shared<Parameter>(
                ParamId{ParameterType::SM, "DB1", LhaID{1}},
                2.0 * a, 0.0, 0.0
            );
            self->store_or_assign(LhaID{1}, out);
        }
    );

    (void)ba.getValue("DB1", C);

    // PARAMS[p] = DB1[c] + 5
    {
        ParamId src_id{ParameterType::SM, "DB1", C};
        ParamId dst_id{ParameterType::SM, "PARAMS", P};

        add_dependent_parameter_local(
            ba,
            dst_id,
            {
                {src_id, ba.getParameter("DB1", C)}
            },
            [src_id](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
                double c = src.get_val(src_id);
                self->set_expected_silent(c + 5.0);
            }
        );
    }

    // DB2[z] = DB1[c] + PARAMS[p]
    add_dependent_block_local(
        ba,
        "DB2",
        {
            {"DB1", ba.at("DB1")},
            {"PARAMS", ba.at("PARAMS")}
        },
        [](const BlockSrc& src, std::shared_ptr<DependentBlock> self) {
            double c = src.get_val("DB1", 1);
            double p = src.get_val("PARAMS", 1);
            auto out = std::make_shared<Parameter>(
                ParamId{ParameterType::SM, "DB2", LhaID{1}},
                c + p, 0.0, 0.0
            );
            self->store_or_assign(LhaID{1}, out);
        }
    );

    assert(std::abs(ba.getValue("BASE", A)   - 10.0) < 1e-12);
    assert(std::abs(ba.getValue("DB1",  C)   - 20.0) < 1e-12);
    assert(std::abs(ba.getValue("PARAMS", P) - 25.0) < 1e-12);
    assert(std::abs(ba.getValue("DB2",  Z)   - 45.0) < 1e-12);

    ba.setValue("BASE", A, 20.0);
    assert(std::abs(ba.getValue("DB1",  C)   - 40.0) < 1e-12);
    assert(std::abs(ba.getValue("PARAMS", P) - 45.0) < 1e-12);
    assert(std::abs(ba.getValue("DB2",  Z)   - 85.0) < 1e-12);

    ba.detach_parameter("PARAMS", P);
    ba.setValue("BASE", A, 7.0);

    assert(std::abs(ba.getValue("DB1",  C)   - 14.0) < 1e-12);
    assert(std::abs(ba.getValue("PARAMS", P) - 45.0) < 1e-12);
    assert(std::abs(ba.getValue("DB2",  Z)   - 59.0) < 1e-12);

    ba.reattach_parameter("PARAMS", P);

    assert(std::abs(ba.getValue("DB1",  C)   - 14.0) < 1e-12);
    assert(std::abs(ba.getValue("PARAMS", P) - 19.0) < 1e-12);
    assert(std::abs(ba.getValue("DB2",  Z)   - 33.0) < 1e-12);

    ba.detach_block("DB1");
    ba.setValue("BASE", A, 100.0);

    assert(std::abs(ba.getValue("DB1",  C)   - 14.0) < 1e-12);
    assert(std::abs(ba.getValue("PARAMS", P) - 19.0) < 1e-12);
    assert(std::abs(ba.getValue("DB2",  Z)   - 33.0) < 1e-12);

    ba.reattach_block("DB1");

    assert(std::abs(ba.getValue("DB1",  C)   - 200.0) < 1e-12);
    assert(std::abs(ba.getValue("PARAMS", P) - 205.0) < 1e-12);
    assert(std::abs(ba.getValue("DB2",  Z)   - 405.0) < 1e-12);
}

int main(){
    std::cout << "== Running BlockAccessor integration tests ==\n";

    auto src1 = std::make_shared<Block>(); src1->blockname = "SRC1";
    auto src2 = std::make_shared<Block>(); src2->blockname = "SRC2";

    LhaID k(100);
    src1->store(k, mkp("SRC1", 100, 1.0));
    src2->store(k, mkp("SRC2", 100, 2.0));

    auto depSum = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{ {"SRC1", src1}, {"SRC2", src2} },
        [k](const auto& blocks, std::shared_ptr<DependentBlock> self){
            double v = blocks.get_val("SRC1",k)
                     + blocks.get_val("SRC2",k);
            if (!self->contains(k)) self->store(k, mkp("SUM", 100, 0.0));
            self->assign(k, v);
        }
    );
    depSum->blockname = "SUM";
    depSum->init();

    auto depPost = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{ {"SUM", depSum} },
        [k](const auto& blocks, std::shared_ptr<DependentBlock> self){
            double base = blocks.get_val("SUM",k);
            if (!self->contains(k)) self->store(k, mkp("POST", 100, 0.0));
            self->assign(k, 3.0 * base);
        }
    );
    depPost->blockname = "POST";
    depPost->init();

    auto acc = std::make_shared<BlockAccessor>();
    acc->emplace(src1->blockname, src1);
    acc->emplace(src2->blockname, src2);
    acc->emplace(depSum->blockname, depSum);
    acc->emplace(depPost->blockname, depPost);

    std::cout << acc->getValue("SUM",  k) << std::endl;
    std::cout << acc->getValue("SRC1",  100) << std::endl;
    std::cout << acc->getValue("SRC2",  100) << std::endl;
    acc->setValue("SRC1", k, 1.0);
    assert(std::abs(acc->getValue("SUM",  k) - 3.0) < 1e-12);
    assert(std::abs(acc->getValue("POST", k) - 9.0) < 1e-12);

    std::cout << acc->getValue("SUM",  k) << std::endl;
    std::cout << acc->getValue("SRC1",  100) << std::endl;
    std::cout << acc->getValue("SRC2",  100) << std::endl;

    acc->setValue("SRC2", k, 5.0);

    std::cout << acc->getValue("SUM",  k) << std::endl;
    std::cout << acc->getValue("SRC1",  100) << std::endl;
    std::cout << acc->getValue("SRC2",  100) << std::endl;
    
    assert(std::abs(acc->getValue("SUM",  k) - 6.0) < 1e-12);
    assert(std::abs(acc->getValue("POST", k) - 18.0) < 1e-12);

    
    acc->at("SUM")->freeze();
    acc->setValue("SRC1", k, 10.0);
    assert(std::abs(acc->getValue("SUM",  k) - 6.0)  < 1e-12);
    assert(std::abs(acc->getValue("POST", k) - 18.0) < 1e-12);

    acc->at("SUM")->unfreeze();
    assert(std::abs(acc->getValue("SUM",  k) - 15.0) < 1e-12);
    assert(std::abs(acc->getValue("POST", k) - 45.0) < 1e-12);

    auto sub = (*acc)[ std::unordered_set<BlockName>{"SRC1", "SUM"} ];
    assert(sub->contains("SRC1") && sub->contains("SUM"));
    assert(!sub->contains("SRC2"));
    assert(!sub->contains("POST"));

    auto extra = std::make_shared<Block>(); extra->blockname = "EXTRA";
    extra->store(k, mkp("EXTRA", 100, 7.0));
    auto acc2 = std::make_shared<BlockAccessor>();
    acc2->emplace("EXTRA", extra);

    auto accU = (acc + acc2);
    assert(accU->contains("EXTRA"));
    assert(std::abs(accU->getValue("EXTRA", k) - 7.0) < 1e-12);
    assert(std::abs(accU->getValue("SUM",   k) - 15.0) < 1e-12);
    assert(std::abs(accU->getValue("POST",  k) - 45.0) < 1e-12);

    auto acc3 = std::make_shared<BlockAccessor>();
    auto src1_alt = std::make_shared<Block>(); src1_alt->blockname = "SRC1";
    src1_alt->store(k, mkp("SRC1", 100, 1.0));
    acc3->emplace("SRC1", src1_alt);

    auto prio = (acc >> acc3);
    assert(prio->contains("SRC1"));
    assert(std::abs(prio->getValue("SRC1", k) - 10.0) < 1e-12);

    prio->remove_item("SRC1", k);
    assert(!prio->has_param("SRC1", k));

    prio->setValue("SRC1", k, 2.0);
    assert(std::abs(prio->getValue("SRC1", k) - 2.0) < 1e-12);

    src1->set_scale(123.0);
    assert(acc->has_scale("SRC1"));
    assert(std::abs(acc->get_scale("SRC1") - 123.0) < 1e-12);

    test_block_accessor_mixed_chain_detach_reattach();
    
    std::cout << "\n BlockAccessor integration suite passed.\n";
    return 0;
}
