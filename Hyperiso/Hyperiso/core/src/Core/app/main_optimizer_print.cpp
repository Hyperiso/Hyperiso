#include <iostream>
#include <chrono>
#include <memory>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <cmath>
#include <iomanip>

#include "Block.h"
#include "BlockAccessor.h"
#include "DependentParameter.h"
#include "ParamOptimizer.h"
#include "SourcesView.h"

static void print_header(const std::string& title) {
    std::cout << "\n==== " << title << " ====\n";
}
static void print_kv(const std::string& k, const std::string& v) {
    std::cout << "  - " << std::left << std::setw(18) << k << " : " << v << "\n";
}
static std::string fmt(double x) {
    std::ostringstream oss; oss << std::setprecision(10) << x; return oss.str();
}
static std::shared_ptr<Parameter> make_param(const std::string& block, const LhaID& id, double v) {
    return std::make_shared<Parameter>(ParamId(block, id), v, 0., 0.);
}
static void dump_block_subset(const std::string& name,
                              const std::shared_ptr<Block>& b,
                              std::initializer_list<int> idxs) {
    std::cout << name << " { ";
    bool first = true;
    for (int i : idxs) {
        if (!first) std::cout << " | ";
        first = false;
        LhaID id(i);
        if (b->contains(id)) std::cout << i << ":" << fmt(b->retrieve(id)->get_val());
        else                 std::cout << i << ":(NA)";
    }
    std::cout << " }\n";
}

static void ensure_ready(const std::shared_ptr<Block>& b,
                         const std::vector<LhaID>& ids,
                         const std::string& tag,
                         int max_tries = 4)
{
    for (int t = 0; t < max_tries; ++t) {
        bool ok = true;
        for (auto& id : ids) {
            if (!b->contains(id)) { ok = false; break; }
        }
        if (ok) {
            std::cout << "[ensure_ready] " << tag << " OK (try " << t << ")\n";
            return;
        }
        if (auto dep = std::dynamic_pointer_cast<DependentBlock>(b)) dep->update();
        else break;
    }
    std::cerr << "[ensure_ready] " << tag << " KO — IDs manquants : ";
    for (auto& id : ids) if (!b->contains(id)) std::cerr << id << " ";
    std::cerr << "\n";
}

int main() {
    using clock = std::chrono::steady_clock;

    const int N   = 12;   
    const int OPS = 2000; 

    const int i1_raw = 5,   i2_raw = 77,  j1_raw = 3,  j2_raw = 10;
    const int i1 = (N ? (i1_raw % N) : 0);
    const int i2 = (N ? (i2_raw % N) : 0);
    const int j1 = (N ? (j1_raw % N) : 0);
    const int j2 = (N ? (j2_raw % N) : 0);

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

    int upd_FUSED_AB = 0, upd_FUSED_BC = 0, upd_SCALED_A = 0, upd_MIXED = 0, upd_REDUCED = 0, upd_META = 0;
    int upd_DP_AZ_1 = 0, upd_DP_AZ_2 = 0, upd_DP_MS_1 = 0, upd_DP_MS_2 = 0, upd_DP_META = 0;

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

    auto AVG = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_A", SRC_A}, {"SRC_B", SRC_B}, {"SRC_C", SRC_C}
        },
        [N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            for (int i = 0; i < N; ++i) {
                LhaID id(i);
                double a = srcs.get_val("SRC_A", id);
                double b = srcs.get_val("SRC_B", id);
                double c = srcs.get_val("SRC_C", id);
                double avg = (a + b + c) / 3.0;
                if (self->contains(id)) self->assign(id, avg);
                else self->store(id, std::make_shared<Parameter>(ParamId(self->get_name(), id), avg, 0., 0.));
            }
        }
    );
    AVG->blockname = "AVG"; AVG->set_scale(1.0); AVG->init();

    auto MIXED = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"SRC_A", SRC_A}, {"FUSED_AB", FUSED_AB}, {"FUSED_BC", FUSED_BC}, {"AVG", AVG}
        },
        [&upd_MIXED, N](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_MIXED;
            for (int i = 0; i < N; ++i) {
                LhaID id(i);
                double a  = srcs.get_val("SRC_A", id);
                double z  = srcs.get_val("FUSED_AB", id);
                double y  = srcs.get_val("FUSED_BC", id);
                double g  = srcs.get_val("AVG", id);
                double m  = 0.5*a + z + 0.1*y - g;
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
    AVG->update();

    MIXED->unfreeze();   MIXED->update();
    REDUCED->unfreeze(); REDUCED->update();

    ensure_ready(FUSED_AB, {LhaID(i1), LhaID(i2)}, "FUSED_AB needs i1/i2");
    ensure_ready(FUSED_BC, {LhaID(i1), LhaID(i2)}, "FUSED_BC needs i1/i2");
    ensure_ready(SCALED_A, {LhaID(j1), LhaID(j2)}, "SCALED_A needs j1/j2");
    ensure_ready(MIXED,    {LhaID(j1), LhaID(j2)}, "MIXED needs j1/j2");

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

    auto DP_AZ_1 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED1", LhaID(1001)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {A_i1, SRC_A->retrieve(LhaID(i1))}, {Z_i1, FUSED_AB->retrieve(LhaID(i1))}
        },
        [A_i1, Z_i1, &upd_DP_AZ_1](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            ++upd_DP_AZ_1; self->set_expected(src.get_val(A_i1) + src.get_val(Z_i1));
        }
    ); DP_AZ_1->init(); DERIVED1->store(LhaID(1001), DP_AZ_1);

    auto DP_AZ_2 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED1", LhaID(1002)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {A_i2, SRC_A->retrieve(LhaID(i2))}, {Z_i2, FUSED_AB->retrieve(LhaID(i2))}
        },
        [A_i2, Z_i2, &upd_DP_AZ_2](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            ++upd_DP_AZ_2; self->set_expected(src.get_val(A_i2) - src.get_val(Z_i2));
        }
    ); DP_AZ_2->init(); DERIVED1->store(LhaID(1002), DP_AZ_2);

    auto DP_MS_1 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED2", LhaID(2001)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {M_j1, MIXED->retrieve(LhaID(j1))}, {S_j1, SCALED_A->retrieve(LhaID(j1))}
        },
        [M_j1, S_j1, &upd_DP_MS_1](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            ++upd_DP_MS_1; self->set_expected(src.get_val(M_j1) + src.get_val(S_j1));
        }
    ); DP_MS_1->init(); DERIVED2->store(LhaID(2001), DP_MS_1);

    auto DP_MS_2 = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED2", LhaID(2002)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            {M_j2, MIXED->retrieve(LhaID(j2))}, {S_j2, SCALED_A->retrieve(LhaID(j2))}
        },
        [M_j2, S_j2, &upd_DP_MS_2](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            ++upd_DP_MS_2; self->set_expected(src.get_val(M_j2) - 0.3*src.get_val(S_j2));
        }
    ); DP_MS_2->init(); DERIVED2->store(LhaID(2002), DP_MS_2);

    auto META = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{
            {"DERIVED1", DERIVED1}, {"REDUCED", REDUCED}
        },
        [&upd_META](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            ++upd_META;
            double v = srcs.get_val("DERIVED1", LhaID(1001))
                     + srcs.get_val("DERIVED1", LhaID(1002))
                     + srcs.get_val("REDUCED",  LhaID(0));
            if (self->contains(LhaID(0))) self->assign(LhaID(0), v);
            else self->store(LhaID(0), std::make_shared<Parameter>(ParamId(self->get_name(), LhaID(0)), v, 0., 0.));
        }
    );
    META->blockname = "META"; META->set_scale(1.0); META->init();
    META->update();

    auto DERIVED_META = std::make_shared<Block>(); DERIVED_META->blockname = "DERIVED_META"; DERIVED_META->set_scale(1.0);
    auto META0_id = ParamId("META", LhaID(0));
    auto S_i1     = SCALED_A->retrieve(LhaID(i1))->get_id();

    auto DP_META = std::make_shared<DependentParameter>(
        ParamId{ParameterType::SM, "DERIVED_META", LhaID(3001)},
        std::unordered_map<ParamId, std::shared_ptr<Parameter>>{
            { META0_id, META->retrieve(LhaID(0)) },
            { S_i1,     SCALED_A->retrieve(LhaID(i1)) }
        },
        [META0_id, S_i1, &upd_DP_META](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            ++upd_DP_META;
            self->set_expected(src.get_val(META0_id) + src.get_val(S_i1));
        }
    ); DP_META->init(); DERIVED_META->store(LhaID(3001), DP_META);


    print_header("INIT (après priming)");
    dump_block_subset("SRC_A",     SRC_A,     {0,1,2,i1,i2});
    dump_block_subset("SRC_B",     SRC_B,     {0,1,2,i1,i2});
    dump_block_subset("SRC_C",     SRC_C,     {0,1,2,i1,i2});
    dump_block_subset("FUSED_AB",  FUSED_AB,  {0,1,2,i1,i2});
    dump_block_subset("FUSED_BC",  FUSED_BC,  {0,1,2,i1,i2});
    dump_block_subset("SCALED_A",  SCALED_A,  {0,1,2,i1,i2});
    dump_block_subset("AVG",       AVG,       {0,1,2,i1,i2});
    dump_block_subset("MIXED",     MIXED,     {0,1,2,j1,j2});
    std::cout << "REDUCED[0] = " << fmt(REDUCED->retrieve(LhaID(0))->get_val()) << "\n";

    auto BA1 = std::make_shared<BlockAccessor>();
    BA1->emplace("SRC_A", SRC_A);
    BA1->emplace("SRC_B", SRC_B);
    BA1->emplace("FUSED_AB", FUSED_AB);
    BA1->emplace("SCALE", SCALE);
    BA1->emplace("SCALED_A", SCALED_A);
    BA1->emplace("AVG", AVG);

    auto BA2 = std::make_shared<BlockAccessor>();
    BA2->emplace("SRC_C", SRC_C);
    BA2->emplace("FUSED_BC", FUSED_BC);
    BA2->emplace("MIXED", MIXED);
    BA2->emplace("REDUCED", REDUCED);
    BA2->emplace("DERIVED1", DERIVED1);
    BA2->emplace("DERIVED2", DERIVED2);
    BA2->emplace("META", META);
    BA2->emplace("DERIVED_META", DERIVED_META);

    ParamOptimizer opt({BA1, BA2});

    auto t0 = clock::now();
    for (int k = 0; k < OPS; ++k) {
        int idx = (N ? (k % N) : 0);
        SRC_A->assign(LhaID(idx), 0.1*k + std::sin(0.001*k));
        SRC_B->assign(LhaID(idx), 0.2*k + std::cos(0.001*k));
        SRC_C->assign(LhaID(idx), 0.05*k);
        if ((k % 200) == 0) SCALE->assign(LhaID(1), 1.25 + 0.0005*k);

        if ((k % (OPS/2)) == 0) {
            print_header(std::string("NAIF checkpoint k=") + std::to_string(k));
            dump_block_subset("FUSED_AB", FUSED_AB, {i1,i2});
            dump_block_subset("FUSED_BC", FUSED_BC, {i1,i2});
            dump_block_subset("SCALED_A", SCALED_A, {i1,i2});
            dump_block_subset("MIXED",    MIXED,    {j1,j2});
            std::cout << "REDUCED[0] = " << fmt(REDUCED->retrieve(LhaID(0))->get_val()) << "\n";
        }
    }
    auto t1 = clock::now();
    auto naive_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    print_header("NAIF — résumé");
    print_kv("upd FUSED_AB",  std::to_string(upd_FUSED_AB));
    print_kv("upd FUSED_BC",  std::to_string(upd_FUSED_BC));
    print_kv("upd SCALED_A",  std::to_string(upd_SCALED_A));
    print_kv("upd MIXED",     std::to_string(upd_MIXED));
    print_kv("upd REDUCED",   std::to_string(upd_REDUCED));
    print_kv("upd DP_AZ_1/2", std::to_string(upd_DP_AZ_1) + "/" + std::to_string(upd_DP_AZ_2));
    print_kv("upd DP_MS_1/2", std::to_string(upd_DP_MS_1) + "/" + std::to_string(upd_DP_MS_2));
    print_kv("upd META",      std::to_string(upd_META));
    print_kv("time (ms)",     std::to_string(naive_ms));
    dump_block_subset("MIXED (post-naif)", MIXED, {j1,j2});
    std::cout << "REDUCED[0] (post-naif) = " << fmt(REDUCED->retrieve(LhaID(0))->get_val()) << "\n";

    upd_FUSED_AB = upd_FUSED_BC = upd_SCALED_A = upd_MIXED = upd_REDUCED = upd_META = 0;
    upd_DP_AZ_1 = upd_DP_AZ_2 = upd_DP_MS_1 = upd_DP_MS_2 = upd_DP_META = 0;

    auto t2 = clock::now();
    for (int k = 0; k < OPS; ++k) {
        int idx = (N ? (k % N) : 0);
        opt.set_value("SRC_A", LhaID(idx), 0.1*k + std::sin(0.002*k));
        opt.set_value("SRC_B", LhaID(idx), 0.2*k + std::cos(0.002*k));
        opt.set_value("SRC_C", LhaID(idx), 0.05*k);
        if ((k % 200) == 0) opt.set_value("SCALE", LhaID(1), 1.30 + 0.0005*k);
    }
    opt.commit(true);
    auto t3 = clock::now();
    auto opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    print_header("OPTIMIZER — résumé");
    print_kv("upd FUSED_AB",  std::to_string(upd_FUSED_AB));
    print_kv("upd FUSED_BC",  std::to_string(upd_FUSED_BC));
    print_kv("upd SCALED_A",  std::to_string(upd_SCALED_A));
    print_kv("upd MIXED",     std::to_string(upd_MIXED));
    print_kv("upd REDUCED",   std::to_string(upd_REDUCED));
    print_kv("upd DP_AZ_1/2", std::to_string(upd_DP_AZ_1) + "/" + std::to_string(upd_DP_AZ_2));
    print_kv("upd DP_MS_1/2", std::to_string(upd_DP_MS_1) + "/" + std::to_string(upd_DP_MS_2));
    print_kv("upd META",      std::to_string(upd_META));
    print_kv("time (ms)",     std::to_string(opt_ms));
    dump_block_subset("MIXED (post-opt)", MIXED, {j1,j2});
    std::cout << "REDUCED[0] (post-opt) = " << fmt(REDUCED->retrieve(LhaID(0))->get_val()) << "\n";

    auto check_point = [&](int raw){
        if (N == 0) return;
        int idx = ((raw % N) + N) % N;
        double a = SRC_A->retrieve(LhaID(idx))->get_val();
        double b = SRC_B->retrieve(LhaID(idx))->get_val();
        double c = SRC_C->retrieve(LhaID(idx))->get_val();
        double z = FUSED_AB->retrieve(LhaID(idx))->get_val();
        double y = FUSED_BC->retrieve(LhaID(idx))->get_val();
        double g = AVG->retrieve(LhaID(idx))->get_val();
        double sf = SCALE->retrieve(LhaID(1))->get_val();
        double s  = SCALED_A->retrieve(LhaID(idx))->get_val();
        double m  = MIXED->retrieve(LhaID(idx))->get_val();
        assert(std::abs(z - (a + 2.0*b)) < 1e-8);
        assert(std::abs(y - (b - c)) < 1e-8);
        assert(std::abs(g - ((a+b+c)/3.0)) < 1e-8);
        assert(std::abs(s - (sf * a)) < 1e-8);
        assert(std::abs(m - (0.5*a + z + 0.1*y - g)) < 1e-8);
    };
    check_point(0); check_point(1); check_point(i1); check_point(i2);

    print_header("FIN — valeurs clés");
    dump_block_subset("SRC_A",     SRC_A,     {i1,i2});
    dump_block_subset("FUSED_AB",  FUSED_AB,  {i1,i2});
    dump_block_subset("SCALED_A",  SCALED_A,  {j1,j2});
    dump_block_subset("MIXED",     MIXED,     {j1,j2});
    std::cout << "REDUCED[0] = " << fmt(REDUCED->retrieve(LhaID(0))->get_val()) << "\n";
    std::cout << "DERIVED1[1001] (A[i1]+Z[i1]) = " << fmt(DERIVED1->retrieve(LhaID(1001))->get_val()) << "\n";
    std::cout << "DERIVED1[1002] (A[i2]-Z[i2]) = " << fmt(DERIVED1->retrieve(LhaID(1002))->get_val()) << "\n";
    std::cout << "DERIVED2[2001] (M[j1]+S[j1]) = " << fmt(DERIVED2->retrieve(LhaID(2001))->get_val()) << "\n";
    std::cout << "DERIVED2[2002] (M[j2]-0.3*S[j2]) = " << fmt(DERIVED2->retrieve(LhaID(2002))->get_val()) << "\n";

    std::cout << "\nOK: maillage complexe — prints + cohérence + bench naïf vs opt.\n";
    return 0;
}
