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

static std::shared_ptr<Parameter> mkp(const std::string& blk, long code, double v){
    return std::make_shared<Parameter>(ParamId{ParameterType::SM, blk, LhaID(code)}, v, 0.0, 0.0);
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

    acc->setValue("SRC1", k, 1.0);
    assert(std::abs(acc->getValue("SUM",  k) - 3.0) < 1e-12);
    assert(std::abs(acc->getValue("POST", k) - 9.0) < 1e-12);

    acc->setValue("SRC2", k, 5.0);
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

    std::cout << "\n BlockAccessor integration suite passed.\n";
    return 0;
}
