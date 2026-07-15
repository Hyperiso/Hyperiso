#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "DependentParameter.h"
#include "Include.h"
#include "SourcesView.h"

static std::shared_ptr<Parameter> mkp(ParameterType type, const std::string& blk, long code, double v){
    return std::make_shared<Parameter>(ParamId{type, blk, LhaID(code)}, v, 0.0, 0.0);
}

static void print_block(BlockAccessor& ba, const std::string& name) {
    std::cout << "\n=== BLOCK " << name << " ===\n";
    if (!ba.contains(name)) {
        std::cout << "(missing)\n";
        return;
    }
    auto blk = ba.at(name);
    auto ids = blk->getAllIDs();
    std::cout << "IDs: ";
    for (auto& id : ids) std::cout << id.to_string() << " ";
    std::cout << "\n";
    for (auto& [id, p] : blk->getItems()) {
        std::cout << "  " << id.to_string() << " = " << p->get_val() << "\n";
    }
}

int main() {
    std::cout << "== Deep dependency graph smoke test ==\n";

    LhaID ID_X(100);
    LhaID ID_P1(200);
    LhaID ID_P2(201);
    LhaID ID_SUM(300);
    LhaID ID_POST(400);

    auto A = std::make_shared<Block>(); A->blockname = "SRC_A";
    auto B = std::make_shared<Block>(); B->blockname = "SRC_B";
    A->store(ID_X, mkp(ParameterType::SM, "SRC_A", 100, 1.0));
    B->store(ID_X, mkp(ParameterType::SM, "SRC_B", 100, 2.0));

    auto acc = std::make_shared<BlockAccessor>();
    acc->emplace("SRC_A", A);
    acc->emplace("SRC_B", B);

    auto depSUM = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_A", A},
            {"SRC_B", B}
        },
        [=](const auto& blocks, std::shared_ptr<DependentBlock> self){
            double v = blocks.get_val("SRC_A", ID_X) + blocks.get_val("SRC_B", ID_X);
            if (!self->contains(ID_SUM)) {
                self->store(ID_SUM, mkp(ParameterType::SM, "SUM", 300, 0.0));
            }
            self->assign(ID_SUM, v);
        }
    );
    depSUM->blockname = "SUM";
    depSUM->init();
    acc->emplace("SUM", depSUM);

    auto DER = std::make_shared<Block>(); DER->blockname = "DERIVED";
    acc->emplace("DERIVED", DER);

    {
        ParamId pid{ParameterType::SM, "DERIVED", ID_P1};
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> srcs;
        // source = param A[ID_X]
        srcs.emplace(ParamId{ParameterType::SM, "SRC_A", ID_X}, A->retrieve(ID_X));

        auto P1 = std::make_shared<DependentParameter>(
            pid, std::move(srcs),
            [=](const ParamSrc& s, std::shared_ptr<DependentParameter> self){
                double a = s.get_val(ParamId{ParameterType::SM, "SRC_A", ID_X});
                self->set_expected(10.0 * a);
            }
        );
        P1->init();
        DER->store_or_assign(ID_P1, P1);
    }

    {
        ParamId pid{ParameterType::SM, "DERIVED", ID_P2};
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> srcs;

        srcs.emplace(ParamId{ParameterType::SM, "SUM", ID_SUM}, depSUM->retrieve(ID_SUM));

        srcs.emplace(ParamId{ParameterType::SM, "DERIVED", ID_P1}, DER->retrieve(ID_P1));

        auto P2 = std::make_shared<DependentParameter>(
            pid, std::move(srcs),
            [=](const ParamSrc& s, std::shared_ptr<DependentParameter> self){
                double sumv = s.get_val(ParamId{ParameterType::SM, "SUM", ID_SUM});
                double p1v  = s.get_val(ParamId{ParameterType::SM, "DERIVED", ID_P1});
                self->set_expected(sumv + p1v);
            }
        );
        P2->init();
        DER->store_or_assign(ID_P2, P2);
    }

    auto depPOST = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SUM", depSUM},
            {"DERIVED", DER}
        },
        [=](const auto& blocks, std::shared_ptr<DependentBlock> self){
            double sumv = blocks.get_val("SUM", ID_SUM);
            double p2v  = blocks.get_val("DERIVED", ID_P2);
            double v = 3.0 * sumv + p2v;
            if (!self->contains(ID_POST)) {
                self->store(ID_POST, mkp(ParameterType::SM, "POST", 400, 0.0));
            }
            self->assign(ID_POST, v);
        }
    );
    depPOST->blockname = "POST";
    depPOST->init();
    acc->emplace("POST", depPOST);


    std::cout << "\n[1] Test lazy via contains() sur blocks dependants\n";
    bool sum_has = acc->at("SUM")->contains(ID_SUM);
    std::cout << "SUM contains(300) = " << sum_has << "\n";
    assert(sum_has && "SUM should materialize on contains()");

    bool post_has = acc->at("POST")->contains(ID_POST);
    std::cout << "POST contains(400) = " << post_has << "\n";
    assert(post_has && "POST should materialize on contains()");

    std::cout << "\n[2] Lecture valeurs (force ensure_up_to_date en cascade)\n";
    double sum0  = acc->getValue("SUM",  ID_SUM);
    double p1_0  = acc->getValue("DERIVED", ID_P1);
    double p2_0  = acc->getValue("DERIVED", ID_P2);
    double post0 = acc->getValue("POST", ID_POST);

    std::cout << "SUM  = " << sum0  << " (expected 1+2=3)\n";
    std::cout << "P1   = " << p1_0  << " (expected 10*A=10)\n";
    std::cout << "P2   = " << p2_0  << " (expected SUM+P1=13)\n";
    std::cout << "POST = " << post0 << " (expected 3*SUM+P2=22)\n";

    assert(std::abs(sum0  - 3.0)  < 1e-12);
    assert(std::abs(p1_0  - 10.0) < 1e-12);
    assert(std::abs(p2_0  - 13.0) < 1e-12);
    assert(std::abs(post0 - 22.0) < 1e-12);

    std::cout << "\n[3] Mutation source: setValue(A=7) => cascade\n";
    acc->setValue("SRC_A", ID_X, 7.0);
    double sum1  = acc->getValue("SUM", ID_SUM);        // 7+2=9
    double p1_1  = acc->getValue("DERIVED", ID_P1);     // 70
    double p2_1  = acc->getValue("DERIVED", ID_P2);     // 79
    double post1 = acc->getValue("POST", ID_POST);      // 3*9+79=106

    std::cout << "SUM  = " << sum1  << " (expected 9)\n";
    std::cout << "P1   = " << p1_1  << " (expected 70)\n";
    std::cout << "P2   = " << p2_1  << " (expected 79)\n";
    std::cout << "POST = " << post1 << " (expected 106)\n";

    assert(std::abs(sum1  - 9.0)   < 1e-12);
    assert(std::abs(p1_1  - 70.0)  < 1e-12);
    assert(std::abs(p2_1  - 79.0)  < 1e-12);
    assert(std::abs(post1 - 106.0) < 1e-12);

    std::cout << "\n[4] Freeze: freeze POST, mutate B, POST doit rester stable\n";
    acc->at("POST")->freeze();
    acc->setValue("SRC_B", ID_X, 10.0);

    double sum2  = acc->getValue("SUM", ID_SUM);
    double post2 = acc->getValue("POST", ID_POST);

    std::cout << "SUM  = " << sum2  << " (expected 17)\n";
    std::cout << "POST = " << post2 << " (expected STILL 106 because frozen)\n";
    assert(std::abs(sum2 - 17.0) < 1e-12);
    assert(std::abs(post2 - 106.0) < 1e-12);

    std::cout << "\n[5] Unfreeze: POST doit se mettre à jour au prochain get\n";
    acc->at("POST")->unfreeze();
    double post3 = acc->getValue("POST", ID_POST);

    std::cout << "POST = " << post3 << " (expected 138)\n";
    assert(std::abs(post3 - 138.0) < 1e-12);

    std::cout << "\n[6] Sanity: getAllIDs/getItems prints\n";
    print_block(*acc, "SRC_A");
    print_block(*acc, "SRC_B");
    print_block(*acc, "SUM");
    print_block(*acc, "DERIVED");
    print_block(*acc, "POST");

    std::cout << "\n[7] Remove source parameter A[100] then check behavior\n";

    acc->remove_item("SRC_A", ID_X);

    bool threw = !acc->at("SUM")->contains(300);
    assert(threw);


    threw = !acc->at("DERIVED")->contains(ID_P1);

    assert(threw);
    std::cout << "\n Deep dependency graph smoke test passed.\n";
    return 0;
}
