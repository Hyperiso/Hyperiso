#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <limits>
#include <chrono>

#include "Block.h"
#include "DependentParameter.h"
#include "SourcesView.h"

/*
================================================================================
Stress-test + micro-bench for deep/wide dependency graphs: Blocks, DependentBlocks,
and DependentParameters, mixed in arbitrary orders/depths.

Graph shape
-----------
WIDTH base blocks at depth 0:
    B0_0 ... B0_{WIDTH-1}

DEPTH levels of dependent blocks:
    D(d,j) for d=1..DEPTH, j=0..WIDTH-1

Block dependencies (DependentBlock)
-----------------------------------
Each DependentBlock D(d,j) depends on TWO blocks from previous level:
    A = level[d-1][j]
    B = level[d-1][(j+1) mod WIDTH]

Important: we intentionally introduce dependencies *through DependentParameters* too.

Base-level DependentParameters
------------------------------
Each base block B0_j also stores one DependentParameter at LhaID(9000 + j):
    BASE_DP(j) = B0_j[1] + 2 * B0_j[2]
This creates:
    param -> DependentParameter (inside a base Block)

DependentBlock recomputation rule (math)
----------------------------------------
For each parameter index i in [1..NPARAMS], define:

    depA = A[9000 + idx(A)]   (this is a DependentParameter stored in A if d-1==0,
                              or a DependentParameter stored in a DependentBlock if d-1>0)
    depB = B[9000 + idx(B)]

    D(d,j)[i] = A[i] + B[i] + 0.1*i + 0.001*(depA + depB)

This creates explicit chains:
    (base params) -> BASE_DP -> DependentBlock params
and also:
    DependentParameter (in source block) -> DependentBlock

DependentParameters stored in DependentBlocks
---------------------------------------------
Each DependentBlock D(d,j) stores two DependentParameters:

  (1) MIX_DP(d,j) at LhaID(1000 + j):
      MIX_DP = A[1] + B[2] + D(d,j)[3]
      (depends on upstream params + a local parameter of the DependentBlock)

  (2) CHAIN_DP(d,j) at LhaID(9000 + j):
      CHAIN_DP = A[9000+idx(A)] + B[9000+idx(B)] + D(d,j)[1]
      (depends on DependentParameters in source blocks + local DependentBlock param)

Thus we have chains like:
    base param -> base DependentParameter -> DependentBlock param -> DependentParameter -> DependentBlock param ...
and many mixed combinations (param/block/depparam/depblock) at arbitrary depths.

What we test
------------
Correctness:
  - We implement a pure C++ expected-value model mirroring the above rules, and we assert
    actual == expected for many random nodes after each update pattern.

Update propagation:
  - Single-leaf updates, batch updates before reads, weird read orders.

Graph operations:
  - freeze/unfreeze semantics (cached until unfreeze).
  - DependentParameter::rebind correctness.

Performance signal (micro-bench):
  - We time loops of (update leaf -> read root) to get a stable notion of update+propagation cost.
  - We also time repeated reads (no updates) to estimate cached access cost.

Notes
-----
This file assumes Include.h provides:
  - enum class ParameterType
  - struct ParamId { ParameterType type; std::string block; int code; ... }
  - struct LhaID with LhaID(int) ctor
Adjust constructors if your Include.h differs.
================================================================================
*/

static constexpr double EPS = 1e-10;

static void assert_near(double got, double exp, const char* what) {
    if (!std::isfinite(got) || !std::isfinite(exp)) {
        std::cerr << "ASSERT_NEAR failed for " << what << " non-finite\n";
        std::abort();
    }

    const double abs_tol = 1e-9;      // plus permissif
    const double rel_tol = 1e-12;     // serré en relatif
    const double diff = std::abs(got - exp);
    const double scale = std::max(1.0, std::abs(exp));

    if (diff > abs_tol && diff > rel_tol * scale) {
        std::cerr << "ASSERT_NEAR failed for " << what
                  << " got=" << got << " expected=" << exp
                  << " diff=" << (got - exp) << "\n";
        std::abort();
    }
}

static std::shared_ptr<DependentParameter> make_dep_param(
    const ParamId& id,
    const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& sources,
    std::function<void(const ParamSrc&, std::shared_ptr<DependentParameter>)> recalc)
{
    auto dp = std::make_shared<DependentParameter>(id, sources, recalc);
    dp->init();
    return dp;
}

static std::shared_ptr<Block> make_base_block(const std::string& name, int n_params, double base)
{
    auto b = std::make_shared<Block>();
    b->blockname = name;
    b->bind_self(b);

    for (int i = 1; i <= n_params; ++i) {
        ParamId pid{ParameterType::SM, name, i};
        auto p = std::make_shared<Parameter>(pid, base + i, /*stat*/0.0, /*syst*/0.0);
        b->store(LhaID(i), p);
    }
    return b;
}

static void add_base_dependent_parameter(std::shared_ptr<Block> b, int base_index_j)
{
    // BASE_DP(j) at LhaID(9000+j) = B0_j[1] + 2*B0_j[2]
    const int code = 9000 + base_index_j;
    ParamId myid{ParameterType::SM, b->blockname, code};

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> psrc;
    psrc[ParamId{ParameterType::SM, b->blockname, 1}] = b->retrieve(LhaID(1));
    psrc[ParamId{ParameterType::SM, b->blockname, 2}] = b->retrieve(LhaID(2));

    auto recalc = [](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
        double v1 = 0.0, v2 = 0.0;
        // raw() gives ParamId->Parameter*, but easiest is sum with weights manually:
        for (auto& [pid, p] : src.raw()) {
            if (!p) continue;
            if (pid.code == LhaID(1)) v1 = p->get_val();
            if (pid.code == LhaID(2)) v2 = p->get_val();
        }
        self->set_expected_silent(v1 + 2.0 * v2);
    };

    auto dp = make_dep_param(myid, psrc, recalc);
    b->store(LhaID(code), dp);
}

static std::shared_ptr<DependentBlock> make_dep_block(
    const std::string& name,
    std::unordered_map<std::string, std::shared_ptr<Block>> sources,
    int n_params,
    int width)
{
    // Rule:
    //   depA = A[9000 + idx(A)]   depB = B[9000 + idx(B)]
    //   self[i] = A[i] + B[i] + 0.1*i + 0.001*(depA + depB)
    //
    // We encode idx(A) from its blockname suffix after "B0_" or "D{d}_".
    // For simplicity, we rely on the fact we stored dep params at 9000+j for every block with column j.
    auto idx_from_name = [](const std::string& bn) -> int {
        // bn looks like "B0_3" or "D7_3" => take part after last '_'
        auto pos = bn.find_last_of('_');
        if (pos == std::string::npos) return 0;
        return std::stoi(bn.substr(pos + 1));
    };

    auto recalc = [n_params, width, idx_from_name](const BlockSrc& src, std::shared_ptr<DependentBlock> self) {
        // identify our two source names in stable manner
        std::vector<std::string> names;
        names.reserve(self->get_source_blocks().size());
        for (auto& [bn, _] : self->get_source_blocks()) names.push_back(bn);
        if (names.size() < 2) return;

        // For robustness: use the first two (in practice we pass exactly 2).
        const std::string& Aname = names[0];
        const std::string& Bname = names[1];
        const int Aj = idx_from_name(Aname) % width;
        const int Bj = idx_from_name(Bname) % width;

        const int depA_id = 9000 + Aj;
        const int depB_id = 9000 + Bj;

        const double depA = src.get_val(Aname, depA_id);
        const double depB = src.get_val(Bname, depB_id);

        for (int i = 1; i <= n_params; ++i) {
            const double a = src.get_val(Aname, i);
            const double b = src.get_val(Bname, i);
            const double val = a + b + 0.1 * i + 0.001 * (depA + depB);

            LhaID id(i);
            if (!self->contains(id)) {
                ParamId pid{ParameterType::SM, self->blockname, i};
                auto p = std::make_shared<Parameter>(pid, 0.0, 0.0, 0.0);
                self->store(id, p);
            }
            self->retrieve(id)->set_expected_silent(val);
        }
    };

    auto db = std::make_shared<DependentBlock>(sources, recalc);
    db->blockname = name;
    db->bind_self(db);
    db->init();
    return db;
}

static void add_depblock_dependent_parameters(
    std::shared_ptr<DependentBlock> db,
    std::shared_ptr<Block> A,
    std::shared_ptr<Block> B,
    int j_col)
{
    // (1) MIX_DP at 1000+j: A[1] + B[2] + D[3]
    {
        ParamId myid{ParameterType::SM, db->blockname, 1000 + j_col};
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> psrc;
        psrc[ParamId{ParameterType::SM, A->blockname, 1}] = A->retrieve(LhaID(1));
        psrc[ParamId{ParameterType::SM, B->blockname, 2}] = B->retrieve(LhaID(2));
        (void)db->retrieve(LhaID(3));
        psrc[ParamId{ParameterType::SM, db->blockname, 3}] = db->retrieve(LhaID(3));

        auto recalc = [](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            double s = 0.0;
            for (auto& [_, p] : src.raw()) if (p) s += p->get_val();
            self->set_expected_silent(s);
        };

        auto dp = make_dep_param(myid, psrc, recalc);
        db->store(LhaID(1000 + j_col), dp);
    }

    // (2) CHAIN_DP at 9000+j: A[9000+Aj] + B[9000+Bj] + D[1]
    // This creates: depparam(source blocks) -> depparam(this block) -> depblock downstream
    {
        ParamId myid{ParameterType::SM, db->blockname, 9000 + j_col};
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> psrc;

        // A and B are previous-level blocks, each must have dep param at 9000+col
        // We assume column index is the suffix after '_' in the block name.
        auto idx_from_name = [](const std::string& bn) -> int {
            auto pos = bn.find_last_of('_');
            if (pos == std::string::npos) return 0;
            return std::stoi(bn.substr(pos + 1));
        };
        const int Aj = idx_from_name(A->blockname);
        const int Bj = idx_from_name(B->blockname);

        psrc[ParamId{ParameterType::SM, A->blockname, 9000 + Aj}] = A->retrieve(LhaID(9000 + Aj));
        psrc[ParamId{ParameterType::SM, B->blockname, 9000 + Bj}] = B->retrieve(LhaID(9000 + Bj));
        (void)db->retrieve(LhaID(1));
        psrc[ParamId{ParameterType::SM, db->blockname, 1}] = db->retrieve(LhaID(1));

        auto recalc = [](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
            double s = 0.0;
            for (auto& [_, p] : src.raw()) if (p) s += p->get_val();
            self->set_expected_silent(s);
        };

        auto dp = make_dep_param(myid, psrc, recalc);
        db->store(LhaID(9000 + j_col), dp);
    }
}

// ----------------------------------------------------------------------------
// Expected-value model mirroring the graph rules
// ----------------------------------------------------------------------------

struct LeafOverride {
    std::unordered_map<long long, double> v;
    static long long key(int j, int i) {
        return (static_cast<long long>(j) << 32) ^ static_cast<unsigned long long>(i);
    }
    void set(int j, int i, double val) { v[key(j,i)] = val; }
    double get(int j, int i, double def) const {
        auto it = v.find(key(j,i));
        return (it == v.end()) ? def : it->second;
    }
};

static double expected_base_val(int base_j, int i, const LeafOverride& ovr) {
    double def = 100.0 * base_j + i;
    return ovr.get(base_j, i, def);
}

static double expected_base_dp(int base_j, const LeafOverride& ovr) {
    // BASE_DP(j) = B0_j[1] + 2*B0_j[2]
    return expected_base_val(base_j, 1, ovr) + 2.0 * expected_base_val(base_j, 2, ovr);
}

static double expected_chain_dp(int d, int j, int width, int nparams, const LeafOverride& ovr); // fwd

static double expected_block_val(int d, int j, int i, int width, int nparams, const LeafOverride& ovr) {
    if (d == 0) return expected_base_val(j, i, ovr);

    int ja = j;
    int jb = (j + 1) % width;

    const double a_i = expected_block_val(d-1, ja, i, width, nparams, ovr);
    const double b_i = expected_block_val(d-1, jb, i, width, nparams, ovr);

    // depA and depB come from chain_dp at previous level blocks
    const double depA = expected_chain_dp(d-1, ja, width, nparams, ovr);
    const double depB = expected_chain_dp(d-1, jb, width, nparams, ovr);

    return a_i + b_i + 0.1 * i + 0.001 * (depA + depB);
}

static double expected_mix_dp(int d, int j, int width, int nparams, const LeafOverride& ovr) {
    // MIX_DP(d,j) = A[1] + B[2] + D(d,j)[3]
    int ja = j;
    int jb = (j + 1) % width;
    const double a1 = expected_block_val(d-1, ja, 1, width, nparams, ovr);
    const double b2 = expected_block_val(d-1, jb, 2, width, nparams, ovr);
    const double self3 = expected_block_val(d, j, 3, width, nparams, ovr);
    return a1 + b2 + self3;
}

static double expected_chain_dp(int d, int j, int width, int nparams, const LeafOverride& ovr) {
    // CHAIN_DP at 9000+j:
    //   if d==0: BASE_DP(j)
    //   else:    CHAIN_DP(d,j) = CHAIN_DP(d-1,j) + CHAIN_DP(d-1,(j+1)) + D(d,j)[1]
    if (d == 0) return expected_base_dp(j, ovr);
    int ja = j;
    int jb = (j + 1) % width;
    const double depA = expected_chain_dp(d-1, ja, width, nparams, ovr);
    const double depB = expected_chain_dp(d-1, jb, width, nparams, ovr);
    const double self1 = expected_block_val(d, j, 1, width, nparams, ovr);
    return depA + depB + self1;
}

// ----------------------------------------------------------------------------
// Graph
// ----------------------------------------------------------------------------

struct Graph {
    int width = 0;
    int depth = 0;
    int nparams = 0;
    std::vector<std::vector<std::shared_ptr<Block>>> levels; // levels[d][j]
};

static Graph build_graph(int width, int depth, int nparams)
{
    Graph g;
    g.width = width;
    g.depth = depth;
    g.nparams = nparams;

    g.levels.emplace_back();
    for (int j = 0; j < width; ++j) {
        auto b = make_base_block("B0_" + std::to_string(j), nparams, 100.0 * j);
        add_base_dependent_parameter(b, j); // adds dep param at 9000+j
        g.levels[0].push_back(b);
    }

    for (int d = 1; d <= depth; ++d) {
        g.levels.emplace_back();
        for (int j = 0; j < width; ++j) {
            auto A = g.levels[d-1][j];
            auto B = g.levels[d-1][(j + 1) % width];

            std::unordered_map<std::string, std::shared_ptr<Block>> src;
            src[A->blockname] = A;
            src[B->blockname] = B;

            auto db = make_dep_block("D" + std::to_string(d) + "_" + std::to_string(j), src, nparams, width);

            // Add MIX_DP and CHAIN_DP inside this DependentBlock.
            add_depblock_dependent_parameters(db, A, B, j);

            g.levels[d].push_back(db);
        }
    }
    return g;
}

// ----------------------------------------------------------------------------
// Test helpers
// ----------------------------------------------------------------------------

static void check_many_values(const Graph& g, const LeafOverride& ovr, std::mt19937_64& rng)
{
    std::uniform_int_distribution<int> ddist(0, g.depth);
    std::uniform_int_distribution<int> jdist(0, g.width - 1);
    std::uniform_int_distribution<int> idist(1, g.nparams);

    // Check random normal params
    for (int t = 0; t < 400; ++t) {
        int d = ddist(rng), j = jdist(rng), i = idist(rng);
        auto b = g.levels[d][j];

        double got = b->retrieve(LhaID(i))->get_val();
        double exp = expected_block_val(d, j, i, g.width, g.nparams, ovr);

        std::string what = "block[" + std::to_string(d) + "][" + std::to_string(j) + "].param[" + std::to_string(i) + "]";
        assert_near(got, exp, what.c_str());
    }

    // Check MIX_DP at random dependent blocks
    std::uniform_int_distribution<int> ddist_dep(1, g.depth);
    for (int t = 0; t < 250; ++t) {
        int d = ddist_dep(rng), j = jdist(rng);
        auto b = g.levels[d][j];

        double got = b->retrieve(LhaID(1000 + j))->get_val();
        double exp = expected_mix_dp(d, j, g.width, g.nparams, ovr);

        std::string what = "block[" + std::to_string(d) + "][" + std::to_string(j) + "].MIX_DP";
        assert_near(got, exp, what.c_str());
    }

    // Check CHAIN_DP at random blocks (including base)
    for (int t = 0; t < 250; ++t) {
        int d = ddist(rng), j = jdist(rng);
        auto b = g.levels[d][j];

        double got = b->retrieve(LhaID(9000 + j))->get_val();
        double exp = expected_chain_dp(d, j, g.width, g.nparams, ovr);

        std::string what = "block[" + std::to_string(d) + "][" + std::to_string(j) + "].CHAIN_DP";
        assert_near(got, exp, what.c_str());
    }
}

static void touch_in_weird_order(const Graph& g)
{
    // Touch nodes in reverse depth order, and interleave dp/params.
    for (int d = g.depth; d >= 0; --d) {
        for (int j = 0; j < g.width; ++j) {
            auto b = g.levels[d][j];
            (void)b->retrieve(LhaID(9000 + j))->get_val(); // CHAIN_DP
            (void)b->retrieve(LhaID(1))->get_val();
            (void)b->retrieve(LhaID(3))->get_val();
            if (d > 0) (void)b->retrieve(LhaID(1000 + j))->get_val(); // MIX_DP
        }
    }
}

static void apply_leaf_updates(Graph& g, LeafOverride& ovr, const std::vector<std::tuple<int,int,double>>& ops)
{
    for (auto& [j,i,v] : ops) {
        g.levels[0][j]->assign(LhaID(i), v);
        ovr.set(j, i, v);
    }
}

static void test_freeze_unfreeze(Graph& g, LeafOverride& ovr)
{
    auto root = g.levels[g.depth][0];
    double before = root->retrieve(LhaID(1))->get_val();

    root->freeze();
    g.levels[0][0]->assign(LhaID(1), 424242.0);
    ovr.set(0, 1, 424242.0);

    double still = root->retrieve(LhaID(1))->get_val();
    assert_near(still, before, "freeze: root[1] cached");

    root->unfreeze();
    double after = root->retrieve(LhaID(1))->get_val();
    double exp = expected_block_val(g.depth, 0, 1, g.width, g.nparams, ovr);
    assert_near(after, exp, "unfreeze: root[1] recompute");
}

static void test_rebind_dependent_parameter(Graph& g, LeafOverride& ovr)
{
    // Pick a DependentParameter (MIX_DP) at deepest level and rebind it:
    // new DP = sourceA[4] + sourceB[5] + sourceA[CHAIN_DP]  (mixing param + depparam)
    int d = g.depth;
    int j = 1 % g.width;

    auto db = std::dynamic_pointer_cast<DependentBlock>(g.levels[d][j]);
    if (!db) { std::cerr << "ERROR: expected DependentBlock\n"; std::abort(); }

    auto dp = std::dynamic_pointer_cast<DependentParameter>(db->retrieve(LhaID(1000 + j)));
    if (!dp) { std::cerr << "ERROR: expected DependentParameter\n"; std::abort(); }

    auto A = g.levels[d-1][j];
    auto B = g.levels[d-1][(j+1)%g.width];

    std::unordered_map<ParamId, std::shared_ptr<Parameter>> new_sources;
    new_sources[ParamId{ParameterType::SM, A->blockname, 4}] = A->retrieve(LhaID(4));
    new_sources[ParamId{ParameterType::SM, B->blockname, 5}] = B->retrieve(LhaID(5));
    new_sources[ParamId{ParameterType::SM, A->blockname, 9000 + j}] = A->retrieve(LhaID(9000 + j));

    auto new_lambda = [](const ParamSrc& src, std::shared_ptr<DependentParameter> self) {
        double s = 0.0;
        for (auto& [_, p] : src.raw()) if (p) s += p->get_val();
        self->set_expected_silent(s);
    };

    double oldv = dp->get_val();
    dp->rebind(new_sources, new_lambda);

    // Mutate an upstream leaf that affects A[4] to force change
    g.levels[0][2]->assign(LhaID(4), 123456.0);
    ovr.set(2, 4, 123456.0);

    const double expA4 = expected_block_val(d-1, j, 4, g.width, g.nparams, ovr);
    const double expB5 = expected_block_val(d-1, (j+1)%g.width, 5, g.width, g.nparams, ovr);
    const double expAchain = expected_chain_dp(d-1, j, g.width, g.nparams, ovr);
    const double exp = expA4 + expB5 + expAchain;

    double newv = dp->get_val();
    if (std::abs(newv - oldv) < 1e-12) {
        std::cerr << "ERROR: rebind dp did not change\n";
        std::abort();
    }
    assert_near(newv, exp, "rebind dp expected");
}

static void microbench(Graph& g, LeafOverride& ovr, std::mt19937_64& rng, bool verbose)
{
    using clock = std::chrono::high_resolution_clock;
    using ns = std::chrono::nanoseconds;

    auto root = g.levels[g.depth][0];

    // Warmup
    for (int k = 0; k < 200; ++k) {
        (void)root->retrieve(LhaID(1))->get_val();
        (void)root->retrieve(LhaID(9000 + 0))->get_val();
        (void)root->retrieve(LhaID(1000 + 0))->get_val();
    }

    std::uniform_int_distribution<int> jdist(0, g.width - 1);
    std::uniform_int_distribution<int> idist(1, g.nparams);
    std::uniform_real_distribution<double> vdist(-1e4, 1e4);

    const int N = 1000;

    // Benchmark A: update leaf then read root (forces propagation).
    auto t0 = clock::now();
    double sink = 0.0;
    for (int t = 0; t < N; ++t) {
        int bj = jdist(rng);
        int pi = idist(rng);
        double nv = vdist(rng);

        g.levels[0][bj]->assign(LhaID(pi), nv);
        ovr.set(bj, pi, nv);

        // Force propagation by reading a few root nodes
        sink += root->retrieve(LhaID(1))->get_val();
        sink += root->retrieve(LhaID(9000 + 0))->get_val();   // CHAIN_DP at root col 0
        sink += root->retrieve(LhaID(1000 + 0))->get_val();   // MIX_DP at root col 0

        if (verbose && (t % 200 == 0)) {
            std::cout << "[bench] t=" << t << " updated B0_" << bj << "[" << pi << "]=" << nv
                      << " -> root[1]=" << root->retrieve(LhaID(1))->get_val() << "\n";
        }
    }
    auto t1 = clock::now();
    auto dt_update_read = std::chrono::duration_cast<ns>(t1 - t0).count();

    // Benchmark B: repeated reads only (cache path).
    auto r0 = clock::now();
    for (int t = 0; t < N * 5; ++t) {
        sink += root->retrieve(LhaID(1))->get_val();
        sink += root->retrieve(LhaID(9000 + 0))->get_val();
        sink += root->retrieve(LhaID(1000 + 0))->get_val();
    }
    auto r1 = clock::now();
    auto dt_reads = std::chrono::duration_cast<ns>(r1 - r0).count();

    std::cout << "\n==== Microbench (rough signal) ====\n";
    std::cout << "Pattern A: (update leaf -> read root[1]+root[CHAIN_DP]+root[MIX_DP]) x " << N << "\n";
    std::cout << "  total: " << (dt_update_read / 1e6) << " ms"
              << " | avg: " << (dt_update_read / (double)N) << " ns/iter\n";
    std::cout << "Pattern B: (read root trio) x " << (N*5) << " (no updates)\n";
    std::cout << "  total: " << (dt_reads / 1e6) << " ms"
              << " | avg: " << (dt_reads / (double)(N*5)) << " ns/iter\n";
    std::cout << "sink=" << sink << " (ignore)\n\n";
}

int main()
{
    constexpr bool VERBOSE = true;

    constexpr int WIDTH = 6;
    constexpr int DEPTH = 8;
    constexpr int NPARAMS = 30;

    std::mt19937_64 rng(123456789ULL);

    std::cout << "Building graph (WIDTH=" << WIDTH << ", DEPTH=" << DEPTH << ", NPARAMS=" << NPARAMS << ")...\n";
    Graph g = build_graph(WIDTH, DEPTH, NPARAMS);
    LeafOverride ovr;

    // Show a small "diagram" sample in output
    std::cout << "Dependency pattern:\n";
    std::cout << "  D(d,j) depends on A=level[d-1][j], B=level[d-1][(j+1)%WIDTH]\n";
    std::cout << "  D(d,j)[i] = A[i] + B[i] + 0.1*i + 0.001*(A[CHAIN_DP] + B[CHAIN_DP])\n";
    std::cout << "  BASE CHAIN_DP in B0_j: B0_j[1] + 2*B0_j[2]\n";
    std::cout << "  MIX_DP in each D(d,j): A[1] + B[2] + D(d,j)[3]\n";
    std::cout << "  CHAIN_DP in each D(d,j): A[CHAIN_DP] + B[CHAIN_DP] + D(d,j)[1]\n\n";

    // 1) Initial correctness sampling
    std::cout << "Test 1: initial correctness (random sampling)...\n";
    check_many_values(g, ovr, rng);

    // 2) Weird read order (lazy / any order)
    std::cout << "Test 2: weird read order then correctness...\n";
    touch_in_weird_order(g);
    check_many_values(g, ovr, rng);

    // 3) Single-leaf update & targeted checks
    std::cout << "Test 3: single-leaf update then checks...\n";
    apply_leaf_updates(g, ovr, { {2, 1, 9999.0} });

    auto root = g.levels[g.depth][0];
    double got_root_1 = root->retrieve(LhaID(1))->get_val();
    double exp_root_1 = expected_block_val(g.depth, 0, 1, g.width, g.nparams, ovr);
    assert_near(got_root_1, exp_root_1, "target: root[1] after leaf update");

    double got_root_chain = root->retrieve(LhaID(9000 + 0))->get_val();
    double exp_root_chain = expected_chain_dp(g.depth, 0, g.width, g.nparams, ovr);
    assert_near(got_root_chain, exp_root_chain, "target: root CHAIN_DP after leaf update");

    double got_root_mix = root->retrieve(LhaID(1000 + 0))->get_val();
    double exp_root_mix = expected_mix_dp(g.depth, 0, g.width, g.nparams, ovr);
    assert_near(got_root_mix, exp_root_mix, "target: root MIX_DP after leaf update");

    check_many_values(g, ovr, rng);

    // 4) Batch updates before reads
    std::cout << "Test 4: batch updates before reads...\n";
    {
        std::vector<std::tuple<int,int,double>> ops;
        std::uniform_int_distribution<int> jdist(0, g.width - 1);
        std::uniform_int_distribution<int> idist(1, g.nparams);
        std::uniform_real_distribution<double> vdist(-1e4, 1e4);

        for (int t = 0; t < 40; ++t) ops.emplace_back(jdist(rng), idist(rng), vdist(rng));
        apply_leaf_updates(g, ovr, ops);

        // Read mid-level dep param first, then root param, to catch ordering bugs.
        auto mid = g.levels[g.depth/2][3];
        double got_mid_chain = mid->retrieve(LhaID(9000 + 3))->get_val();
        double exp_mid_chain = expected_chain_dp(g.depth/2, 3, g.width, g.nparams, ovr);
        assert_near(got_mid_chain, exp_mid_chain, "batch: mid CHAIN_DP");

        double got_root_3 = root->retrieve(LhaID(3))->get_val();
        double exp_root_3 = expected_block_val(g.depth, 0, 3, g.width, g.nparams, ovr);
        assert_near(got_root_3, exp_root_3, "batch: root[3]");

        check_many_values(g, ovr, rng);
    }

    // 5) Freeze/unfreeze
    std::cout << "Test 5: freeze/unfreeze...\n";
    test_freeze_unfreeze(g, ovr);

    // 6) Rebind dependent parameter
    std::cout << "Test 6: DependentParameter::rebind...\n";
    test_rebind_dependent_parameter(g, ovr);

    // 7) Microbench timing output
    std::cout << "Bench: timing signal...\n";
    microbench(g, ovr, rng, VERBOSE);

    std::cout << "ALL TESTS PASSED ✅\n";
    return 0;
}