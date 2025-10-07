#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>

#include "ParameterRouter.h"
#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"
#include "General.h"

static std::shared_ptr<Parameter> mkp(const BlockName& blk, long code, double v){
    return std::make_shared<Parameter>(ParamId{ParameterType::SM, blk, LhaID(code)}, v, 0.0, 0.0);
}

int main(){
    std::cout << "== Running ParameterRouter integration tests ==\n";

    auto acc = std::make_shared<BlockAccessor>();

    auto b_mass  = std::make_shared<Block>();  b_mass->blockname  = "MASS";
    auto b_gauge = std::make_shared<Block>();  b_gauge->blockname = "GAUGE";
    auto b_alpha = std::make_shared<Block>();  b_alpha->blockname = "ALPHA"; // BSM

    b_mass->store(LhaID(25), mkp("MASS", 25, 0.1));
    b_mass->store(LhaID(35), mkp("MASS", 35, 0.2));

    b_gauge->store(LhaID(2),  mkp("GAUGE", 2,  1.0));
    b_gauge->store(LhaID(99), mkp("GAUGE", 99, 2.0));

    b_alpha->store(LhaID(0), mkp("ALPHA", 0, 3.14));

    acc->emplace("MASS",  b_mass);
    acc->emplace("GAUGE", b_gauge);
    acc->emplace("ALPHA", b_alpha);


    {
        auto p25 = acc->at("MASS")->retrieve(LhaID(25));
        auto p35 = acc->at("MASS")->retrieve(LhaID(35));

        auto t25 = ParamRouter::GetType("MASS", LhaID(25));
        auto t35 = ParamRouter::GetType("MASS", LhaID(35));

        p25->set_owner(t25);
        p35->set_owner(t35);

        assert(p25->get_id().type == ParameterType::SM);
        assert(p35->get_id().type == ParameterType::BSM);
    }

    {
        auto p2  = acc->at("GAUGE")->retrieve(LhaID(2));
        auto p99 = acc->at("GAUGE")->retrieve(LhaID(99));

        p2->set_owner(ParamRouter::GetType("GAUGE", LhaID(2)));
        p99->set_owner(ParamRouter::GetType("GAUGE", LhaID(99)));

        assert(p2->get_id().type  == ParameterType::SM);
        assert(p99->get_id().type == ParameterType::BSM);
    }

    {
        auto pa = acc->at("ALPHA")->retrieve(LhaID(0));
        pa->set_owner(ParamRouter::GetType("ALPHA", LhaID(0)));
        assert(pa->get_id().type == ParameterType::BSM);
    }

    {
        auto v_mass = ParamRouter::GetType("MASS");
        std::unordered_set<ParameterType> s(v_mass.begin(), v_mass.end());
        assert(s.count(ParameterType::SM)  == 1);
        assert(s.count(ParameterType::BSM) == 1);

        auto v_alpha = ParamRouter::GetType("ALPHA");
        std::unordered_set<ParameterType> sa(v_alpha.begin(), v_alpha.end());
        assert(sa.size() == 1 && sa.count(ParameterType::BSM) == 1);
    }

    {
        auto sm_blocks = ParamRouter::GetOwnedBlocks(ParameterType::SM);
        auto bsm_blocks = ParamRouter::GetOwnedBlocks(ParameterType::BSM);
        assert(sm_blocks.count("MASS") == 1 && sm_blocks.count("GAUGE") == 1);
        assert(bsm_blocks.count("MASS") == 1 && bsm_blocks.count("ALPHA") == 1);
    }

    {
        std::vector<BlockName> present = {"MASS", "GAUGE", "ALPHA"};
        auto customs = ParameterBlockRepartition::filter_custom_blocks(present);
        std::unordered_set<BlockName> S(customs.begin(), customs.end());

        assert(S.count("MASS") == 1);
        assert(S.count("GAUGE") == 1);

        assert(S.count("ALPHA") == 0);
    }

    std::cout << "\n ParameterRouter integration suite passed.\n";
    return 0;
}
