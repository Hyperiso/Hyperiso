#include <iostream>
#include <chrono>
#include <memory>
#include <unordered_map>
#include <cassert>

#include "Block.h"
#include "BlockAccessor.h"
#include "DependentParameter.h"
#include "ParamOptimizer.h"
#include "SourcesView.h"

static std::shared_ptr<Parameter> make_param(const std::string& block, const LhaID& id, double v) {
    return std::make_shared<Parameter>(ParamId(block, id), v, 0., 0.);
}

int main() {
    using clock = std::chrono::steady_clock;

    const int N   = 400;
    const int OPS = 20000;
    const int id_sum = 3;

    auto SRC_A = std::make_shared<Block>(); SRC_A->blockname = "SRC_A"; SRC_A->set_scale(1.0);
    auto SRC_B = std::make_shared<Block>(); SRC_B->blockname = "SRC_B"; SRC_B->set_scale(1.0);
    for (int i = 0; i < N; ++i) {
        LhaID id(i);
        SRC_A->store(id, make_param("SRC_A", id, 1.0*i));
        SRC_B->store(id, make_param("SRC_B", id, 2.0*i));
    }

    int fused_updates = 0;
    int sum_updates   = 0;

    auto make_fused = [&](auto& fused_ptr){
        fused_ptr = std::make_shared<DependentBlock>(
            std::unordered_map<std::string, std::shared_ptr<Block>>{
                {"SRC_A", SRC_A},
                {"SRC_B", SRC_B}
            },
            [&fused_updates, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
                ++fused_updates;
                for (int i = 0; i < N; ++i) {
                    LhaID id(i);
                    double a = srcs.get_val("SRC_A", id);
                    double b = srcs.get_val("SRC_B", id);
                    double z = a + b;
                    if (self->contains(id)) self->assign(id, z);
                    else self->store(id, std::make_shared<Parameter>(ParamId(self->get_name(), id), z, 0., 0.));
                }
            }
        );
        fused_ptr->blockname = "FUSED"; fused_ptr->set_scale(1.0);
        fused_ptr->init();
        SRC_A->addObserver(fused_ptr);
        SRC_B->addObserver(fused_ptr);
    };

    std::shared_ptr<DependentBlock> FUSED;
    make_fused(FUSED);

    auto DERIVED = std::make_shared<Block>(); DERIVED->blockname = "DERIVED"; DERIVED->set_scale(1.0);
    auto srcA_sum = SRC_A->retrieve(LhaID(id_sum));
    auto srcB_sum = SRC_B->retrieve(LhaID(id_sum));
    const ParamId a_id = srcA_sum->get_id();
    const ParamId b_id = srcB_sum->get_id();

    auto SUM = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED", LhaID(1000)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            { a_id, srcA_sum },
            { b_id, srcB_sum },
        },
        [a_id, b_id, &sum_updates](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            ++sum_updates;
            self->set_expected(src.get_val(a_id) + src.get_val(b_id));
        }
    );
    SUM->init();
    DERIVED->store(LhaID(1000), SUM);

    auto BA1 = std::make_shared<BlockAccessor>();
    BA1->emplace("SRC_A", SRC_A);
    BA1->emplace("FUSED", FUSED);
    auto BA2 = std::make_shared<BlockAccessor>();
    BA2->emplace("SRC_B", SRC_B);
    BA2->emplace("DERIVED", DERIVED);
    ParamOptimizer opt(std::vector<std::shared_ptr<BlockAccessor>>{BA1, BA2});

    //sloooow
    fused_updates = 0; sum_updates = 0;
    auto t0 = clock::now();
    for (int k = 0; k < OPS; ++k) {
        int i = k % N;
        SRC_A->assign(LhaID(i), double(k) * 0.1);
        SRC_B->assign(LhaID(i), double(k) * 0.2);
    }
    auto t1 = clock::now();
    auto naive_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::cout << "[NAIF] updates FUSED = " << fused_updates
              << ", updates SUM = " << sum_updates
              << ", time = " << naive_ms << " ms\n";

    {
        int i = 17;
        double a = SRC_A->retrieve(LhaID(i))->get_val();
        double b = SRC_B->retrieve(LhaID(i))->get_val();
        double z = FUSED->retrieve(LhaID(i))->get_val();
        assert(std::abs(z - (a+b)) < 1e-12);
        assert(std::abs(SUM->get_val() - (srcA_sum->get_val()+srcB_sum->get_val())) < 1e-12);
    }

    //Niiiiiice
    fused_updates = 0; sum_updates = 0;
    auto t2 = clock::now();
    for (int k = 0; k < OPS; ++k) {
        int i = k % N;
        opt.set_value("SRC_A", LhaID(i), double(k) * 0.1);
        opt.set_value("SRC_B", LhaID(i), double(k) * 0.2);
    }
    opt.commit(true);
    auto t3 = clock::now();
    auto opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    std::cout << "[OPT ] updates FUSED = " << fused_updates
              << ", updates SUM = " << sum_updates
              << ", time = " << opt_ms << " ms\n";

    {
        int i = 123;
        double a = SRC_A->retrieve(LhaID(i))->get_val();
        double b = SRC_B->retrieve(LhaID(i))->get_val();
        double z = FUSED->retrieve(LhaID(i))->get_val();
        assert(std::abs(z - (a+b)) < 1e-12);
        assert(std::abs(SUM->get_val() - (srcA_sum->get_val()+srcB_sum->get_val())) < 1e-12);
    }

    // opt, store+remove on new keys ()
    fused_updates = 0; sum_updates = 0;
    for (int j = 0; j < 100; ++j) {
        LhaID new_id(N + j);
        opt.set_value("SRC_A", new_id, 42.0); // store
        opt.remove("SRC_A", new_id);          // sera ignoré (clé inexistante avant commit)
    }
    opt.commit(true);

    {
        int i = 7;
        double a = SRC_A->retrieve(LhaID(i))->get_val();
        double b = SRC_B->retrieve(LhaID(i))->get_val();
        double z = FUSED->retrieve(LhaID(i))->get_val();
        assert(std::abs(z - (a+b)) < 1e-12);
    }

    std::cout << "OK: cohérence et bench terminés.\n";
    return 0;
}
