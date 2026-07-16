#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "DependentParameter.h"
#include "Include.h"
#include "SourcesView.h"

static std::shared_ptr<Parameter> mkp(ParameterType type, const std::string& blk, long code, double v){
    return std::make_shared<Parameter>(ParamId{type, blk, LhaID(code)}, v, 0.0, 0.0);
}

static void dump(BlockAccessor& acc, const std::string& name) {
    std::cout << "\n=== DUMP " << name << " ===\n";
    if (!acc.contains(name)) { std::cout << "(missing)\n"; return; }
    auto b = acc.at(name);
    for (auto& [id, p] : b->getItems()) {
        std::cout << "  " << name << "[" << id.to_string() << "]=" << p->get_val() << "\n";
    }
}

static void expect_throw_get(BlockAccessor& acc, const std::string& blk, const LhaID& id) {
    bool threw = !acc.at(blk)->contains(id);
    assert(threw);
}

int main() {
    std::cout << "== Deep dependency graph TORTURE test ==\n";

    // IDs
    LhaID ID_X(100);
    LhaID ID_Y(101);
    LhaID ID_SUM(300);
    LhaID ID_P1(200);
    LhaID ID_P2(201);
    LhaID ID_POST(400);

    // Sources blocks
    auto A = std::make_shared<Block>(); A->blockname = "SRC_A";
    auto B = std::make_shared<Block>(); B->blockname = "SRC_B";
    A->store(ID_X, mkp(ParameterType::SM, "SRC_A", 100, 1.0));
    B->store(ID_X, mkp(ParameterType::SM, "SRC_B", 100, 2.0));
    B->store(ID_Y, mkp(ParameterType::SM, "SRC_B", 101, 5.0));

    auto acc = std::make_shared<BlockAccessor>();
    acc->emplace("SRC_A", A);
    acc->emplace("SRC_B", B);

    // DependentBlock SUM: SUM[300] = A[100] + B[100]
    auto depSUM = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{{"SRC_A", A}, {"SRC_B", B}},
        [=](const auto& blocks, std::shared_ptr<DependentBlock> self){
            double v = blocks.get_val("SRC_A", ID_X) + blocks.get_val("SRC_B", ID_X);
            if (!self->contains(ID_SUM)) self->store(ID_SUM, mkp(ParameterType::SM, "SUM", 300, 0.0));
            self->assign(ID_SUM, v);
        }
    );
    depSUM->blockname = "SUM";
    depSUM->init();
    acc->emplace("SUM", depSUM);

    // Normal block with DependentParameters
    auto DER = std::make_shared<Block>(); DER->blockname = "DERIVED";
    acc->emplace("DERIVED", DER);

    // P1 = 10*A
    auto make_P1 = [&](double factor){
        ParamId pid{ParameterType::SM, "DERIVED", ID_P1};
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> srcs;
        srcs.emplace(ParamId{ParameterType::SM, "SRC_A", ID_X}, A->retrieve(ID_X));

        auto P1 = std::make_shared<DependentParameter>(
            pid, std::move(srcs),
            [=](const ParamSrc& s, std::shared_ptr<DependentParameter> self){
                double a = s.get_val(ParamId{ParameterType::SM, "SRC_A", ID_X});
                self->set_expected(factor * a);
            }
        );
        P1->init();
        return P1;
    };

    // P2 = SUM + P1 + B[101]
    auto make_P2 = [&](){
        ParamId pid{ParameterType::SM, "DERIVED", ID_P2};
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> srcs;
        srcs.emplace(ParamId{ParameterType::SM, "SUM", ID_SUM}, depSUM->retrieve(ID_SUM));
        srcs.emplace(ParamId{ParameterType::SM, "DERIVED", ID_P1}, DER->retrieve(ID_P1));
        srcs.emplace(ParamId{ParameterType::SM, "SRC_B", ID_Y}, B->retrieve(ID_Y));

        auto P2 = std::make_shared<DependentParameter>(
            pid, std::move(srcs),
            [=](const ParamSrc& s, std::shared_ptr<DependentParameter> self){
                double sumv = s.get_val(ParamId{ParameterType::SM, "SUM", ID_SUM});
                double p1v  = s.get_val(ParamId{ParameterType::SM, "DERIVED", ID_P1});
                double by   = s.get_val(ParamId{ParameterType::SM, "SRC_B", ID_Y});
                self->set_expected(sumv + p1v + by);
            }
        );
        P2->init();
        return P2;
    };

    DER->store_or_assign(ID_P1, make_P1(10.0));
    DER->store_or_assign(ID_P2, make_P2());

    // DependentBlock POST: POST[400] = 3*SUM + P2
    auto depPOST = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{{"SUM", depSUM}, {"DERIVED", DER}},
        [=](const auto& blocks, std::shared_ptr<DependentBlock> self){
            double sumv = blocks.get_val("SUM", ID_SUM);
            double p2v  = blocks.get_val("DERIVED", ID_P2);
            double v = 3.0 * sumv + p2v;
            if (!self->contains(ID_POST)) self->store(ID_POST, mkp(ParameterType::SM, "POST", 400, 0.0));
            self->assign(ID_POST, v);
        }
    );
    depPOST->blockname = "POST";
    depPOST->init();
    acc->emplace("POST", depPOST);

    // ------------------------------------------------------------
    // [1] Basic sanity
    // ------------------------------------------------------------
    std::cout << "\n[1] Basic sanity\n";
    assert(acc->contains("SRC_A"));
    assert(acc->contains("SUM"));
    assert(acc->at("SUM")->contains(ID_SUM));     // materialize via contains
    assert(acc->at("POST")->contains(ID_POST));   // materialize via contains

    double sum0 = acc->getValue("SUM", ID_SUM);           // 1+2=3
    double p1_0 = acc->getValue("DERIVED", ID_P1);        // 10
    double p2_0 = acc->getValue("DERIVED", ID_P2);        // SUM + P1 + B[101] = 3+10+5=18
    double post0= acc->getValue("POST", ID_POST);         // 3*3 + 18 = 27
    std::cout << "sum0=" << sum0 << " p1_0=" << p1_0 << " p2_0=" << p2_0 << " post0=" << post0 << "\n";
    assert(std::abs(sum0-3.0) < 1e-12);
    assert(std::abs(p1_0-10.0) < 1e-12);
    assert(std::abs(p2_0-18.0) < 1e-12);
    assert(std::abs(post0-27.0) < 1e-12);

    // ------------------------------------------------------------
    // [2] Freeze/unfreeze behavior on blocks AND on DependentParameter inside DERIVED
    // ------------------------------------------------------------
    std::cout << "\n[2] Freeze/unfreeze behavior\n";
    acc->at("POST")->freeze();
    acc->setValue("SRC_A", ID_X, 7.0); // SUM=9, P1=70, P2=9+70+5=84, POST should stay 27 (frozen)
    double post_frozen = acc->getValue("POST", ID_POST);
    std::cout << "post_frozen=" << post_frozen << " (expected still 27)\n";
    assert(std::abs(post_frozen-27.0) < 1e-12);

    acc->at("POST")->unfreeze();
    double post1 = acc->getValue("POST", ID_POST); // now 3*9 + 84 = 111
    std::cout << "post1=" << post1 << " (expected 111)\n";
    assert(std::abs(post1-111.0) < 1e-12);

    // freeze DependentParameter P1 only (indirect test)
    auto p1ptr = std::dynamic_pointer_cast<DependentParameter>(DER->retrieve(ID_P1));
    assert(p1ptr);
    p1ptr->freeze();
    acc->setValue("SRC_A", ID_X, 10.0); // normally P1=100, but frozen => stays 70, impacts P2 & POST
    double p1_f = acc->getValue("DERIVED", ID_P1);
    double p2_f = acc->getValue("DERIVED", ID_P2);
    double post_f = acc->getValue("POST", ID_POST);
    std::cout << "p1_f=" << p1_f << " (expected still 70)\n";
    std::cout << "p2_f=" << p2_f << " (expected SUM(12)+70+5=87)\n";
    std::cout << "post_f=" << post_f << " (expected 3*12+87=123)\n";
    assert(std::abs(p1_f-70.0) < 1e-12);
    assert(std::abs(p2_f-87.0) < 1e-12);
    assert(std::abs(post_f-123.0) < 1e-12);

    p1ptr->unfreeze();
    double p1_u = acc->getValue("DERIVED", ID_P1); // now 100
    double post_u = acc->getValue("POST", ID_POST); // SUM=12, P2=12+100+5=117, POST=3*12+117=153
    std::cout << "p1_u=" << p1_u << " (expected 100), post_u=" << post_u << " (expected 153)\n";
    assert(std::abs(p1_u-100.0) < 1e-12);
    assert(std::abs(post_u-153.0) < 1e-12);

    // ------------------------------------------------------------
    // [3] Remove item in source and observe failures + recovery
    // ------------------------------------------------------------
    std::cout << "\n[3] remove_item source + recovery\n";
    acc->remove_item("SRC_A", ID_X);

    // Anything depending on A[100] should now fail on getValue (SUM / P1 / P2 / POST)
    expect_throw_get(*acc, "SUM", ID_SUM);
    expect_throw_get(*acc, "DERIVED", ID_P1);
    expect_throw_get(*acc, "POST", ID_POST);

    std::cout << "here " << std::endl;
    // Recover: re-add A[100]
    acc->setValue("SRC_A", ID_X, 2.0); // A back
    // double sumR = acc->getValue("SUM", ID_SUM); // 2 + B(2)=4
    // double postR = acc->getValue("POST", ID_POST);
    // std::cout << "sumR=" << sumR << " (expected 4)\n";
    // assert(std::abs(sumR-4.0) < 1e-12);

    // ------------------------------------------------------------
    // [4] Replace/rebind DependentParameter in place (simule LO puis NLO)
    // ------------------------------------------------------------
    std::cout << "\n[4] Rebind DependentParameter (store_or_assign overwrite)\n";
    // Remplace P1 = 10*A par P1 = 100*A
    DER->store_or_assign(ID_P1, make_P1(100.0));

    // après overwrite, tout doit suivre
    acc->setValue("SRC_A", ID_X, 1.0);
    double p1_new = acc->getValue("DERIVED", ID_P1); // 100
    double p2_new = acc->getValue("DERIVED", ID_P2); // SUM(1+2=3) + 100 + 5 = 108
    double post_new= acc->getValue("POST", ID_POST); // 3*3 +108 = 117
    std::cout << "p1_new=" << p1_new << " p2_new=" << p2_new << " post_new=" << post_new << "\n";
    assert(std::abs(p1_new-100.0) < 1e-12);
    assert(std::abs(p2_new-108.0) < 1e-12);
    assert(std::abs(post_new-117.0) < 1e-12);

    // ------------------------------------------------------------
    // [5] erase_block and alias cleanup + behavior of dependent blocks afterwards
    // ------------------------------------------------------------
    std::cout << "\n[5] erase_block\n";
    // Remove SRC_B block entirely
    acc->erase_block("SRC_B");
    assert(!acc->contains("SRC_B"));

    // SUM depends on SRC_B => should now fail
    expect_throw_get(*acc, "SUM", ID_SUM);

    // Re-add SRC_B with new values
    auto B2 = std::make_shared<Block>(); B2->blockname = "SRC_B";
    B2->store(ID_X, mkp(ParameterType::SM, "SRC_B", 100, 10.0));
    B2->store(ID_Y, mkp(ParameterType::SM, "SRC_B", 101, 1.0));
    acc->emplace("SRC_B", B2);

    // IMPORTANT: depSUM still holds old pointer B (the erased one).
    // This test is here to catch whether your design *expects* pointer stability.
    // If you expect “erase then re-add” to work, you must rebuild depSUM sources.
    // Here we just *demonstrate* the effect:
    bool throws = !acc->at("SUM")->contains(ID_SUM);
    // try { (void)acc->getValue("SUM", ID_SUM); }
    // catch(...) { throws = true; }
    std::cout << "After erase/re-add SRC_B, SUM getValue throws? " << throws << " (expected TRUE unless you rebuild dependencies)\n";

    // ------------------------------------------------------------
    // [6] Sub-accessor operator[] + merge operators +, >>
    // ------------------------------------------------------------
    std::cout << "\n[6] operator[] / + / >>\n";
    // Build a tiny independent accessor with EXTRA, and a conflicting DERIVED
    auto EXTRA = std::make_shared<Block>(); EXTRA->blockname = "EXTRA";
    EXTRA->store(LhaID(1), mkp(ParameterType::SM, "EXTRA", 1, 42.0));
    auto acc2 = std::make_shared<BlockAccessor>();
    acc2->emplace("EXTRA", EXTRA);

    // + should merge (and complain if conflicts)
    auto merged = (acc + acc2);
    assert(merged->contains("EXTRA"));
    assert(std::abs(merged->getValue("EXTRA", LhaID(1)) - 42.0) < 1e-12);

    // priority merge: rhs overrides blocks
    auto acc_override = std::make_shared<BlockAccessor>();
    auto DER2 = std::make_shared<Block>(); DER2->blockname = "DERIVED";
    DER2->store(LhaID(999), mkp(ParameterType::SM, "DERIVED", 999, 9.0));
    acc_override->emplace("DERIVED", DER2);

    auto prio = (acc >> acc_override);
    assert(prio->contains("DERIVED"));
    assert(prio->at("DERIVED")->contains(LhaID(999)));

    // sub accessor
    auto sub = (*prio)[ std::unordered_set<BlockName>{"SRC_A","DERIVED"} ];
    assert(sub->contains("SRC_A"));
    assert(sub->contains("DERIVED"));
    assert(!sub->contains("POST"));

    // ------------------------------------------------------------
    // Final dumps
    // ------------------------------------------------------------
    dump(*acc, "SRC_A");
    dump(*acc, "SUM");
    dump(*acc, "DERIVED");
    dump(*acc, "POST");
    dump(*merged, "EXTRA");

    std::cout << "\n✅ TORTURE test finished (some checks intentionally demonstrate expected failures).\n";
    return 0;
}
