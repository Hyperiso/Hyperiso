#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_set>

#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "Include.h"

static std::shared_ptr<Parameter> mkp(const std::string& blk, long code, double v){
    return std::make_shared<Parameter>(ParamId{ParameterType::SM, blk, LhaID(code)}, v, 0.0, 0.0);
}

int main() {
    std::cout << "== Running BlockAccessor unit tests ==\n";

    auto bA = std::make_shared<Block>(); bA->blockname = "GAUGE";
    auto bB = std::make_shared<Block>(); bB->blockname = "MASS";

    BlockAccessor ba;
    ba.emplace(bA->blockname, bA);
    ba.emplace(bB->blockname, bB);

    assert(ba.contains("GAUGE"));
    assert(ba.contains("MASS"));
    assert(!ba.contains("YUKAWA"));

    LhaID k1(11), k2(22);
    assert(!ba.has_param("GAUGE", k1));
    assert(!ba.has_param("MASS",  k2));

    ba.setValue("GAUGE", k1, 1.25);
    assert(ba.has_param("GAUGE", k1));
    assert(std::abs(ba.getValue("GAUGE", k1) - 1.25) < 1e-12);

    ba.setValue("GAUGE", k1, 2.75);
    assert(std::abs(ba.getValue("GAUGE", k1) - 2.75) < 1e-12);

    assert(bA->getItems().size() == 1);

    auto pB_k2 = mkp("MASS", 22, 9.0);
    ba.setParameter("MASS", k2, pB_k2);
    assert(ba.has_param("MASS", k2));
    auto got_ptr = ba.getParameter("MASS", k2);

    assert(got_ptr.get() == pB_k2.get());
    assert(std::abs(ba.getValue("MASS", k2) - 9.0) < 1e-12);

    auto valsA = ba.getAllValues("GAUGE");
    assert(valsA.size() == 1);
    assert(valsA.begin()->first == k1);
    assert(std::abs(valsA.begin()->second - 2.75) < 1e-12);

    auto valsB = ba.getAllValues("MASS");
    assert(valsB.size() == 1);
    assert(valsB.begin()->first == k2);
    assert(std::abs(valsB.begin()->second - 9.0) < 1e-12);

    auto names = ba.get_block_names();
    assert(names.size() == 2);
    assert(names.contains("GAUGE"));
    assert(names.contains("MASS"));

    ba.remove_item("MASS", k2);
    assert(!ba.has_param("MASS", k2));

    // ba.remove_item("DOES_NOT_EXIST", LhaID(1));

    {
        auto& ref = ba.at("GAUGE");
        assert(ref.get() == bA.get());

        const BlockAccessor& cba = ba;
        auto& cref = cba.at("GAUGE");
        assert(cref.get() == bA.get());
    }

    assert(!ba.has_scale("GAUGE"));
    bA->set_scale(246.0);
    assert(ba.has_scale("GAUGE"));
    assert(std::abs(ba.get_scale("GAUGE") - 246.0) < 1e-12);

    try {
        (void)ba.getValue("YUKAWA", LhaID(1));
        assert(false && "getValue should have thrown");
    } catch (const std::invalid_argument&) {

    }
    try {
        (void)ba.getParameter("YUKAWA", LhaID(1));
        assert(false && "getParameter should have thrown");
    } catch (const std::invalid_argument&) {

    }

    auto rhs = std::make_shared<BlockAccessor>();
    auto br = std::make_shared<Block>(); br->blockname = "GAUGE";
    rhs->emplace("GAUGE", br);
    LhaID k3(33);
    rhs->setValue("GAUGE", k3, 1.0);

    auto lhs = std::make_shared<BlockAccessor>();
    auto bl = std::make_shared<Block>(); bl->blockname = "GAUGE";
    lhs->emplace("GAUGE", bl);
    lhs->setValue("GAUGE", k3, 123.0);

    auto bx = std::make_shared<Block>(); bx->blockname = "X";
    rhs->emplace("X", bx);
    rhs->setValue("X", LhaID(7), 7.0);

    auto merged = (lhs >> rhs);
    assert(merged->contains("GAUGE") && merged->contains("X"));
    assert(std::abs(merged->getValue("GAUGE", k3) - 123.0) < 1e-12); 
    assert(std::abs(merged->getValue("X", LhaID(7)) - 7.0) < 1e-12);

    auto sub = ba[ std::unordered_set<BlockName>{"GAUGE"} ];
    assert(sub->contains("GAUGE"));
    assert(!sub->contains("MASS"));
    assert(std::abs(sub->getValue("GAUGE", k1) - 2.75) < 1e-12);

    std::cout << "\n BlockAccessor unit suite passed.\n";
    return 0;
}
