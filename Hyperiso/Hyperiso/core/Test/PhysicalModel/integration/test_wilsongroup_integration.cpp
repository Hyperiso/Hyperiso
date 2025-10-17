#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <map>

#include "WilsonGroup.h"
#include "Wilson.h"
#include "BWilsonGroup.h"
#include "IParameterProxy.h"
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "Parameter.h"
#include "Include.h"

class DummyBoolAPI : public ICoreAPI<bool> { public: bool get() override { return false; } };
class DummyStringAPI : public ICoreAPI<std::string> { public: std::string get() override { return "SM"; } };
class DummyPathAPI : public ICoreAPI<fs::path> { public: fs::path get() override { return fs::path("/dev/null"); } };

class SpyProxy : public IParameterProxy<std::string, LhaID> {
public:
    mutable std::string last_block;
    mutable LhaID      last_id{0};
    double ret = 5.0;
    scalar_t operator()(const std::string& b, const LhaID& id) const override { last_block=b; last_id=id; return ret; }
    bool exist(const std::string&, const LhaID&) const override { return true; }
    double get_scale(const std::string&) const override { return 0.0; }
};

class SpyComposer : public IBlockComposer {
public:
    struct R { ParamId target; std::unordered_set<ParamId> sources; };
    std::vector<R> calls;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& id,
                           const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc&) override {
        calls.push_back({id, src});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

int main() {
    std::cout << "== Wilson GROUP INTEGRATION ==\n";

    auto proxy   = std::make_shared<SpyProxy>();
    auto comp    = std::make_shared<SpyComposer>();
    auto use_mty = std::make_shared<DummyBoolAPI>();
    auto mname   = std::make_shared<DummyStringAPI>();
    auto mpath   = std::make_shared<DummyPathAPI>();
    WilsonGroupAdapterConfig cfg(proxy, comp, use_mty, mname, mpath);

    auto grp = std::make_shared<BCoefficientGroup>(cfg, /*force_sm=*/false);
    grp->init(QCDOrder::LO);
    assert(grp->size() == 10);
    assert(comp->calls.size() >= 10);

    const auto match_blk = GroupMapper::str(WGroup::B, ScaleType::MATCHING);
    bool saw_c7_lo = false;
    for (const auto& r : comp->calls) {
        if (r.target == ParamId{match_blk, grp->at("C7")->get_lhaid(QCDOrder::LO)}) {
            ParamId need{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)};
            assert(r.sources.count(need) == 1);
            saw_c7_lo = true;
        }
    }
    assert(saw_c7_lo);

    proxy->ret = 9.0;
    auto vr = grp->get_running_coefficient("C10", "LO", ContributionType::SM, WilsonBasis::B_STANDARD);
    assert(std::abs(vr.real() - 9.0) < 1e-12);
    auto had_blk = GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD);
    assert(proxy->last_block == had_blk);
    assert(proxy->last_id == grp->at("C10")->id(QCDOrder::LO, ContributionType::SM));

    std::cout << "INTEGRATION OK\n";
    return 0;
}
