#include <iostream>
#include <chrono>
#include <memory>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <cmath>

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

    const int N   = 600;
    const int OPS = 300;
    const int i1 = 5, i2 = 77, j1 = 123, j2 = 451;

    auto SRC_A = std::make_shared<Block>(); SRC_A->blockname = "SRC_A"; SRC_A->set_scale(1.0);
    auto SRC_B = std::make_shared<Block>(); SRC_B->blockname = "SRC_B"; SRC_B->set_scale(1.0);
    auto SRC_C = std::make_shared<Block>(); SRC_C->blockname = "SRC_C"; SRC_C->set_scale(1.0);

    auto SCALE = std::make_shared<Block>(); SCALE->blockname = "SCALE"; SCALE->set_scale(1.0);
    SCALE->store(LhaID(1), make_param("SCALE", LhaID(1), 1.25)); // scale factor

    for (int i = 0; i < N; ++i) {
        LhaID id(i);
        SRC_A->store(id, make_param("SRC_A", id, 1.0*i));
        SRC_B->store(id, make_param("SRC_B", id, 2.0*i));
        SRC_C->store(id, make_param("SRC_C", id, 0.5*i));
    }

    int upd_FUSED_AB = 0, upd_FUSED_BC = 0, upd_SCALED_A = 0, upd_MIXED = 0, upd_REDUCED = 0;
    int upd_DP1 = 0, upd_DP2 = 0, upd_DP3 = 0, upd_DP4 = 0;


    auto FUSED_AB = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_A", SRC_A}, {"SRC_B", SRC_B}
        },
        [&upd_FUSED_AB, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_FUSED_AB;
            for (int i = 0; i < N; ++i) {
                LhaID id(i);
                double a = srcs.get_val("SRC_A", id);
                double b = srcs.get_val("SRC_B", id);
                double z = a + 2.0*b;
                if (self->contains(id)) self->assign(id, z);
                else self->store(id, std::make_shared<Parameter>(ParamId(self->get_name(), id), z, 0., 0.));
            }
        }
    );
    FUSED_AB->blockname = "FUSED_AB"; FUSED_AB->set_scale(1.0); FUSED_AB->init();

    auto FUSED_BC = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_B", SRC_B}, {"SRC_C", SRC_C}
        },
        [&upd_FUSED_BC, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_FUSED_BC;
            for (int i = 0; i < N; ++i) {
                LhaID id(i);
                double b = srcs.get_val("SRC_B", id);
                double c = srcs.get_val("SRC_C", id);
                double y = b - c;
                if (self->contains(id)) self->assign(id, y);
                else self->store(id, std::make_shared<Parameter>(ParamId(self->get_name(), id), y, 0., 0.));
            }
        }
    );
    FUSED_BC->blockname = "FUSED_BC"; FUSED_BC->set_scale(1.0); FUSED_BC->init();

    auto SCALED_A = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_A", SRC_A}, {"SCALE", SCALE}
        },
        [&upd_SCALED_A, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_SCALED_A;
            double sf = srcs.get_val("SCALE", 1);
            for (int i = 0; i < N; ++i) {
                LhaID id(i);
                double a = srcs.get_val("SRC_A", id);
                double s = sf * a;
                if (self->contains(id)) self->assign(id, s);
                else self->store(id, std::make_shared<Parameter>(ParamId(self->get_name(), id), s, 0., 0.));
            }
        }
    );
    SCALED_A->blockname = "SCALED_A"; SCALED_A->set_scale(1.0); SCALED_A->init();

    auto MIXED = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_A", SRC_A}, {"FUSED_AB", FUSED_AB}, {"FUSED_BC", FUSED_BC}
        },
        [&upd_MIXED, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_MIXED;
            for (int i = 0; i < N; ++i) {
                LhaID id(i);
                double a  = srcs.get_val("SRC_A", id);
                double z  = srcs.get_val("FUSED_AB", id);
                double y  = srcs.get_val("FUSED_BC", id);
                double m  = 0.5*a + z + 0.1*y;
                if (self->contains(id)) self->assign(id, m);
                else self->store(id, std::make_shared<Parameter>(ParamId(self->get_name(), id), m, 0., 0.));
            }
        }
    );
    MIXED->blockname = "MIXED"; MIXED->set_scale(1.0); MIXED->init();

    auto REDUCED = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"MIXED", MIXED}
        },
        [&upd_REDUCED, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_REDUCED;
            double sum = 0.0;
            for (int i = 0; i < N; ++i) sum += srcs.get_val("MIXED", LhaID(i));
            if (self->contains(LhaID(0))) self->assign(LhaID(0), sum);
            else self->store(LhaID(0), std::make_shared<Parameter>(ParamId(self->get_name(), LhaID(0)), sum, 0., 0.));
        }
    );
    REDUCED->blockname = "REDUCED"; REDUCED->set_scale(1.0); REDUCED->init();


    MIXED->freeze();
    REDUCED->freeze();

    FUSED_AB->update(); 
    FUSED_BC->update(); 
    SCALED_A->update(); 

    MIXED->unfreeze();
    MIXED->update();
    REDUCED->unfreeze();
    REDUCED->update();

    upd_FUSED_AB = upd_FUSED_BC = upd_SCALED_A = upd_MIXED = upd_REDUCED = 0;

    auto DERIVED1 = std::make_shared<Block>(); DERIVED1->blockname = "DERIVED1"; DERIVED1->set_scale(1.0);
    auto DERIVED2 = std::make_shared<Block>(); DERIVED2->blockname = "DERIVED2"; DERIVED2->set_scale(1.0);

    auto A_i1 = SRC_A->retrieve(LhaID(i1))->get_id();
    auto A_i2 = SRC_A->retrieve(LhaID(i2))->get_id();
    auto Z_i1 = FUSED_AB->retrieve(LhaID(i1))->get_id();
    auto Z_i2 = FUSED_AB->retrieve(LhaID(i2))->get_id();
    auto M_j1 = MIXED->retrieve(LhaID(j1))->get_id();
    auto M_j2 = MIXED->retrieve(LhaID(j2))->get_id();
    auto S_j1 = SCALED_A->retrieve(LhaID(j1))->get_id();
    auto S_j2 = SCALED_A->retrieve(LhaID(j2))->get_id();

    auto DP1 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED1", LhaID(1001)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {A_i1, SRC_A->retrieve(LhaID(i1))}, {Z_i1, FUSED_AB->retrieve(LhaID(i1))}
        },
        [A_i1, Z_i1, &upd_DP1](const ParamSrc& src, std::shared_ptr<DependentParameter> self){
            ++upd_DP1; self->set_expected(src.get_val(A_i1) + src.get_val(Z_i1));
        }
    ); DP1->init(); DERIVED1->store(LhaID(1001), DP1);

    auto DP2 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED1", LhaID(1002)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {A_i2, SRC_A->retrieve(LhaID(i2))}, {Z_i2, FUSED_AB->retrieve(LhaID(i2))}
        },
        [A_i2, Z_i2, &upd_DP2](const ParamSrc& src, std::shared_ptr<DependentParameter> self){
            ++upd_DP2; self->set_expected(src.get_val(A_i2) - src.get_val(Z_i2));
        }
    ); DP2->init(); DERIVED1->store(LhaID(1002), DP2);

    auto DP3 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED2", LhaID(2001)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {M_j1, MIXED->retrieve(LhaID(j1))}, {S_j1, SCALED_A->retrieve(LhaID(j1))}
        },
        [M_j1, S_j1, &upd_DP3](const ParamSrc& src, std::shared_ptr<DependentParameter> self){
            ++upd_DP3; self->set_expected(src.get_val(M_j1) + src.get_val(S_j1));
        }
    ); DP3->init(); DERIVED2->store(LhaID(2001), DP3);

    auto DP4 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED2", LhaID(2002)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {M_j2, MIXED->retrieve(LhaID(j2))}, {S_j2, SCALED_A->retrieve(LhaID(j2))}
        },
        [M_j2, S_j2, &upd_DP4](const ParamSrc& src, std::shared_ptr<DependentParameter> self){
            ++upd_DP4; self->set_expected(src.get_val(M_j2) - 0.3*src.get_val(S_j2));
        }
    ); DP4->init(); DERIVED2->store(LhaID(2002), DP4);

    auto BA1 = std::make_shared<BlockAccessor>();
    BA1->emplace("SRC_A", SRC_A);
    BA1->emplace("SRC_B", SRC_B);
    BA1->emplace("FUSED_AB", FUSED_AB);
    BA1->emplace("SCALE", SCALE);
    BA1->emplace("SCALED_A", SCALED_A);

    auto BA2 = std::make_shared<BlockAccessor>();
    BA2->emplace("SRC_C", SRC_C);
    BA2->emplace("FUSED_BC", FUSED_BC);
    BA2->emplace("MIXED", MIXED);
    BA2->emplace("REDUCED", REDUCED);
    BA2->emplace("DERIVED1", DERIVED1);
    BA2->emplace("DERIVED2", DERIVED2);

    ParamOptimizer opt({BA1, BA2});

    upd_FUSED_AB = upd_FUSED_BC = upd_SCALED_A = upd_MIXED = upd_REDUCED = 0;
    upd_DP1 = upd_DP2 = upd_DP3 = upd_DP4 = 0;

    auto t0 = clock::now();
    for (int k = 0; k < OPS; ++k) {
        int i = k % N;
        SRC_A->assign(LhaID(i), 0.1*k + std::sin(0.001*k));
        SRC_B->assign(LhaID(i), 0.2*k + std::cos(0.001*k));
        SRC_C->assign(LhaID(i), 0.05*k);

        if ((k % 1000) == 0) SCALE->assign(LhaID(1), 1.0 + 0.001*k);

        if ((k % 5000) == 0) {
            LhaID nid(N + (k/5000));
            SRC_A->store(nid, make_param("SRC_A", nid, 42.0 + k*1e-3));
            SRC_A->notifyObservers(); // propager l’ajout d’ID
        }
    }
    auto t1 = clock::now();
    auto naive_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::cout << "[NAIF] "
              << "upd(FUSED_AB)=" << upd_FUSED_AB << ", "
              << "upd(FUSED_BC)=" << upd_FUSED_BC << ", "
              << "upd(SCALED_A)=" << upd_SCALED_A << ", "
              << "upd(MIXED)=" << upd_MIXED << ", "
              << "upd(REDUCED)=" << upd_REDUCED << " | "
              << "upd(DP1,DP2,DP3,DP4)=" << upd_DP1 << "," << upd_DP2 << "," << upd_DP3 << "," << upd_DP4
              << " | time=" << naive_ms << " ms\n";

    auto check_point = [&](int idx){
        double a = SRC_A->retrieve(LhaID(idx))->get_val();
        double b = SRC_B->retrieve(LhaID(idx))->get_val();
        double c = SRC_C->retrieve(LhaID(idx))->get_val();
        double z = FUSED_AB->retrieve(LhaID(idx))->get_val();
        double y = FUSED_BC->retrieve(LhaID(idx))->get_val();
        double sf = SCALE->retrieve(LhaID(1))->get_val();
        double s  = SCALED_A->retrieve(LhaID(idx))->get_val();
        double m  = MIXED->retrieve(LhaID(idx))->get_val();
        assert(std::abs(z - (a + 2.0*b)) < 1e-9);
        assert(std::abs(y - (b - c)) < 1e-9);
        assert(std::abs(s - (sf * a)) < 1e-9);
        assert(std::abs(m - (0.5*a + z + 0.1*y)) < 1e-9);
    };
    check_point(3); check_point(17); check_point(123);

    upd_FUSED_AB = upd_FUSED_BC = upd_SCALED_A = upd_MIXED = upd_REDUCED = 0;
    upd_DP1 = upd_DP2 = upd_DP3 = upd_DP4 = 0;

    auto t2 = clock::now();
    for (int k = 0; k < OPS; ++k) {
        int i = k % N;
        opt.set_value("SRC_A", LhaID(i), 0.1*k + std::sin(0.001*k));
        opt.set_value("SRC_B", LhaID(i), 0.2*k + std::cos(0.001*k));
        opt.set_value("SRC_C", LhaID(i), 0.05*k);
        if ((k % 1000) == 0) opt.set_value("SCALE", LhaID(1), 1.0 + 0.001*k);
        if ((k % 5000) == 0) {
            LhaID nid(N + (k/5000));
            opt.set_value("SRC_A", nid, 42.0 + k*1e-3);
        }
    }
    opt.commit(true);
    auto t3 = clock::now();
    auto opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    std::cout << "[OPT ] "
              << "upd(FUSED_AB)=" << upd_FUSED_AB << ", "
              << "upd(FUSED_BC)=" << upd_FUSED_BC << ", "
              << "upd(SCALED_A)=" << upd_SCALED_A << ", "
              << "upd(MIXED)=" << upd_MIXED << ", "
              << "upd(REDUCED)=" << upd_REDUCED << " | "
              << "upd(DP1,DP2,DP3,DP4)=" << upd_DP1 << "," << upd_DP2 << "," << upd_DP3 << "," << upd_DP4
              << " | time=" << opt_ms << " ms\n";

    check_point(7); check_point(111); check_point(482);
    {
        double a1 = SRC_A->retrieve(LhaID(i1))->get_val();
        double z1 = FUSED_AB->retrieve(LhaID(i1))->get_val();
        double a2 = SRC_A->retrieve(LhaID(i2))->get_val();
        double z2 = FUSED_AB->retrieve(LhaID(i2))->get_val();
        double m1 = MIXED->retrieve(LhaID(j1))->get_val();
        double s1 = SCALED_A->retrieve(LhaID(j1))->get_val();
        double m2 = MIXED->retrieve(LhaID(j2))->get_val();
        double s2 = SCALED_A->retrieve(LhaID(j2))->get_val();

        assert(std::abs(DERIVED1->retrieve(LhaID(1001))->get_val() - (a1 + z1)) < 1e-9);
        assert(std::abs(DERIVED1->retrieve(LhaID(1002))->get_val() - (a2 - z2)) < 1e-9);
        assert(std::abs(DERIVED2->retrieve(LhaID(2001))->get_val() - (m1 + s1)) < 1e-9);
        assert(std::abs(DERIVED2->retrieve(LhaID(2002))->get_val() - (m2 - 0.3*s2)) < 1e-9);
    }

    std::cout << "OK: test complexe — cohérence + bench.\n";
    return 0;
}
