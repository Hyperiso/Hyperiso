#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

#include "WilsonManager.h"

#include "CoefficientGroupBuilder.h"
#include "WilsonCoefficientRegistry.h"
#include "GroupMapper.h"
#include "WilsonBlockNames.h"
#include "wcoef_ids.hpp"

#include "WilsonGroupAdapterConfig.h"
#include "IMartyWilsonProxy.h"
#include "InterpretedParam.h"

struct RecordingBlockComposer : IBlockComposer {
    struct ParamCall {
        ParamId dest;
        std::unordered_set<ParamId> sources;
    };
    std::vector<ParamCall> param_calls;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& pid,
                           const std::unordered_set<ParamId>& deps,
                           const DepParamUpdateFunc&) override
    {
        param_calls.push_back({pid, deps});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

template<typename T>
struct DummyCoreAPI : ICoreAPI<T> {
    T value;
    explicit DummyCoreAPI(const T& v) : value(v) {}
    T get() override { return value; }
};

struct DummyScaleSetter : IParamSetter<ScaleType> {
    void set(double) override {}
    void switch_param(ScaleType) override {}
};

static int qcd_index(QCDOrder o) {
    switch (o) {
        case QCDOrder::LO:   return 0;
        case QCDOrder::NLO:  return 1;
        case QCDOrder::NNLO: return 2;
    }
    return 0;
}

static std::string lhaid_key(const LhaID& id) {
    auto p = id.get_parts(); // [custom, base, qcd, part]
    return std::to_string(p[0]) + "_" + std::to_string(p[1]) + "_" +
           std::to_string(p[2]) + "_" + std::to_string(p[3]);
}

static std::unordered_set<std::string> dep_keys(const std::unordered_set<ParamId>& deps) {
    std::unordered_set<std::string> out;
    for (const auto& d : deps) {
        out.insert(lhaid_key(d.code));
    }
    return out;
}

static bool has_call(const RecordingBlockComposer& rec,
                     const ParamId& dest,
                     const std::unordered_set<ParamId>& sources)
{
    const auto want = dep_keys(sources);
    for (const auto& c : rec.param_calls) {
        if (lhaid_key(c.dest.code) != lhaid_key(dest.code)) continue;
        if (c.dest.block != dest.block) continue;
        if (c.dest.type  != dest.type)  continue;
        if (dep_keys(c.sources) == want) return true;
    }
    return false;
}

struct FWExistProxy : IParameterProxy<std::string, LhaID> {
    std::unordered_map<std::string, std::unordered_set<std::string>> exist_keys;

    scalar_t operator()(const std::string&, const LhaID&) const override {
        return 4.0 * PI;
    }

    bool exist(const std::string& block, const LhaID& id) const override {
        auto it = exist_keys.find(block);
        if (it == exist_keys.end()) return false;
        return it->second.contains(lhaid_key(id));
    }

    double get_scale(const std::string&) const override {
        return 1.0;
    }
};

static WilsonGroupAdapterConfig make_adapters(
    const std::shared_ptr<IParameterProxy<std::string, LhaID>>& wilson_proxy,
    const std::shared_ptr<IBlockComposer>& iblock_c,
    const std::shared_ptr<ICoreAPI<bool>>& use_marty)
{
    auto marty_model_name = std::make_shared<DummyCoreAPI<std::string>>(std::string{""});
    auto marty_model_path = std::make_shared<DummyCoreAPI<fs::path>>(fs::path{});
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;

    return WilsonGroupAdapterConfig{
        wilson_proxy,
        iblock_c,
        use_marty,
        marty_model_name,
        marty_model_path,
        marty_proxy
    };
}

struct Harness {
    CoefficientRegistry reg;
    CoefficientGroupBuilder builder;

    std::shared_ptr<RecordingBlockComposer> rec;
    std::shared_ptr<FWExistProxy>           fwproxy;
    std::shared_ptr<DummyCoreAPI<bool>>     use_marty;
    std::shared_ptr<DummyCoreAPI<bool>>     has_input;
    std::shared_ptr<DummyCoreAPI<Model>>    model_api;
    std::shared_ptr<DummyScaleSetter>       scale_setter;

    PortsConfig ports;
    CoefficientManager mgr;

    WGroupId gid;
    std::string groupName;
    std::string fw_block;
    std::string sm_block;
    std::string final_block;

    explicit Harness(Model model)
        : reg()
        , builder(reg)
        , rec(std::make_shared<RecordingBlockComposer>())
        , fwproxy(std::make_shared<FWExistProxy>())
        , use_marty(std::make_shared<DummyCoreAPI<bool>>(false))
        , has_input(std::make_shared<DummyCoreAPI<bool>>(true))
        , model_api(std::make_shared<DummyCoreAPI<Model>>(model))
        , scale_setter(std::make_shared<DummyScaleSetter>())
        , ports(rec, fwproxy, use_marty, has_input, model_api, scale_setter)
        , mgr(ports)
        , gid(GroupMapper::to_id(WGroup::B))
        , groupName(GroupMapper::str(gid))
        , fw_block(WilsonBlockNames::fwcoef())
        , sm_block(WilsonBlockNames::sm_matching(gid))
        , final_block(WilsonBlockNames::matching(gid))
    {
        register_B(reg);

        {
            auto adapters = make_adapters(fwproxy, rec, use_marty);
            BuildContext ctx{adapters, Model::SM, Backend::Builtin, ContributionType::SM, gid};
            auto grp = builder.build(ctx);
            assert(grp && "B group must be buildable");
            mgr.registerCoefficientGroup(groupName, grp);
        }

        ports.build_group = [&](WGroupId gid2,
                                Model model2,
                                bool /*marty_backend*/,
                                ContributionType ctype,
                                std::string storage_block) -> std::shared_ptr<CoefficientGroup>
        {
            auto adapters = make_adapters(fwproxy, rec, use_marty);
            BuildContext ctx{adapters, model2, Backend::Builtin, ctype, gid2};
            auto g = builder.build(ctx);
            assert(g && "build_group must succeed in tests");
            g->set_matching_storage_block(storage_block);
            return g;
        };
    }

    void reset() {
        rec->param_calls.clear();
        fwproxy->exist_keys.clear();
    }

    ParamId pid(const std::string& block, const std::string& coefName, QCDOrder o, int part) {
        auto w = WCoefMapper::enum_of(WCoefMapper::enum_elt(coefName)).value();
        auto base = WCoefMapper::flha_base(w);
        LhaID id(base.first, base.second, qcd_index(o), part);
        return ParamId{ParameterType::WILSON, block, id};
    }
};

static void test_tot_only() {
    Harness H(Model::SUSY);
    H.reset();

    H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C7", QCDOrder::LO, 2).code));

    H.mgr.init_group_matching(H.groupName, "LO");

    auto fw_tot    = H.pid(H.fw_block,    "C7", QCDOrder::LO, 2);
    auto final_sm  = H.pid(H.final_block, "C7", QCDOrder::LO, 0);
    auto final_bsm = H.pid(H.final_block, "C7", QCDOrder::LO, 1);
    auto final_tot = H.pid(H.final_block, "C7", QCDOrder::LO, 2);

    assert(has_call(*H.rec, final_tot, {fw_tot}));

    assert(has_call(*H.rec, final_bsm, {final_tot, final_sm}));
}

static void test_sm_only_rule_bsm0_tot_eq_sm() {
    Harness H(Model::SUSY);
    H.reset();

    H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C7", QCDOrder::LO, 0).code));

    H.mgr.init_group_matching(H.groupName, "LO");

    auto fw_sm     = H.pid(H.fw_block,    "C7", QCDOrder::LO, 0);
    auto sm_inter  = H.pid(H.sm_block,    "C7", QCDOrder::LO, 0);

    auto final_sm  = H.pid(H.final_block, "C7", QCDOrder::LO, 0);
    auto final_bsm = H.pid(H.final_block, "C7", QCDOrder::LO, 1);
    auto final_tot = H.pid(H.final_block, "C7", QCDOrder::LO, 2);

    assert(has_call(*H.rec, sm_inter, {fw_sm}));

    assert(has_call(*H.rec, final_bsm, {}));

    assert(has_call(*H.rec, final_tot, {final_sm}));
}

static void test_tot_and_bsm_no_sm() {
    Harness H(Model::SUSY);
    H.reset();

    H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C7", QCDOrder::LO, 2).code));
    H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C7", QCDOrder::LO, 1).code));

    H.mgr.init_group_matching(H.groupName, "LO");

    auto fw_bsm   = H.pid(H.fw_block,   "C7", QCDOrder::LO, 1);
    auto fw_tot   = H.pid(H.fw_block,   "C7", QCDOrder::LO, 2);
    auto sm_inter = H.pid(H.sm_block,   "C7", QCDOrder::LO, 0);

    auto final_bsm = H.pid(H.final_block, "C7", QCDOrder::LO, 1);
    auto final_tot = H.pid(H.final_block, "C7", QCDOrder::LO, 2);

    assert(has_call(*H.rec, sm_inter, {fw_tot, fw_bsm}));

    assert(has_call(*H.rec, final_bsm, {fw_bsm}));
    assert(has_call(*H.rec, final_tot, {fw_tot}));
}

static void test_bsm_only_tot_eq_sm_plus_bsm() {
    Harness H(Model::SUSY);
    H.reset();

    H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C7", QCDOrder::LO, 1).code));

    H.mgr.init_group_matching(H.groupName, "LO");

    auto fw_bsm    = H.pid(H.fw_block,    "C7", QCDOrder::LO, 1);
    auto final_sm  = H.pid(H.final_block, "C7", QCDOrder::LO, 0);
    auto final_bsm = H.pid(H.final_block, "C7", QCDOrder::LO, 1);
    auto final_tot = H.pid(H.final_block, "C7", QCDOrder::LO, 2);

    assert(has_call(*H.rec, final_bsm, {fw_bsm}));

    assert(has_call(*H.rec, final_tot, {final_sm, final_bsm}));
}

static void test_realistic_fwcoef_block_subset() {
    Harness H(Model::SUSY);
    H.reset();

    for (auto o : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
        H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C7", o, 2).code));
    }
    for (auto o : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
        H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C8", o, 2).code));
    }
    for (auto o : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
        H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C1", o, 2).code));
    }
    for (auto o : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
        H.fwproxy->exist_keys[H.fw_block].insert(lhaid_key(H.pid(H.fw_block, "C2", o, 2).code));
    }

    H.mgr.init_group_matching(H.groupName, "LO");

    {
        auto fw_tot    = H.pid(H.fw_block,    "C7", QCDOrder::LO, 2);
        auto final_sm  = H.pid(H.final_block, "C7", QCDOrder::LO, 0);
        auto final_bsm = H.pid(H.final_block, "C7", QCDOrder::LO, 1);
        auto final_tot = H.pid(H.final_block, "C7", QCDOrder::LO, 2);

        assert(has_call(*H.rec, final_tot, {fw_tot}));
        assert(has_call(*H.rec, final_bsm, {final_tot, final_sm}));
    }

    {
        auto final_sm  = H.pid(H.final_block, "C3", QCDOrder::LO, 0);
        auto final_bsm = H.pid(H.final_block, "C3", QCDOrder::LO, 1);
        auto final_tot = H.pid(H.final_block, "C3", QCDOrder::LO, 2);

        assert(has_call(*H.rec, final_bsm, {}));
        assert(has_call(*H.rec, final_tot, {final_sm}));
    }
}

int main() {
    std::cout << "== has_wilson_input + FWCOEF wiring tests ==\n";

    test_tot_only();
    test_sm_only_rule_bsm0_tot_eq_sm();
    test_tot_and_bsm_no_sm();
    test_bsm_only_tot_eq_sm_plus_bsm();
    test_realistic_fwcoef_block_subset();

    std::cout << "All FWCOEF wiring tests passed.\n";
    return 0;
}
