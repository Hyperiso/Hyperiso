#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "WilsonGroup.h"
#include "Wilson.h"
#include "BWilson.h"
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
    double ret = 123.0;

    scalar_t operator()(const std::string& b, const LhaID& id) const override {
        last_block = b; last_id = id; return ret;
    }
    bool exist(const std::string&, const LhaID&) const override { return true; }
    double get_scale(const std::string&) const override { return 0.0; }
};

struct ParamRecord { ParamId target; std::unordered_set<ParamId> sources; };

class SpyComposer : public IBlockComposer {
public:
    std::vector<ParamRecord> params;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& id,
                           const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc&) override {
        params.push_back({id, src});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

class TestGroupNoInit : public CoefficientGroup {
public:
    TestGroupNoInit(WGroup gid, WilsonGroupAdapterConfig cfg, ContributionType t = ContributionType::SM)
        : CoefficientGroup(cfg)
    {
        this->id = GroupMapper::to_id(gid);
        this->wilson_type = t;
        this->block_name = GroupMapper::str(gid, ScaleType::MATCHING);
    }

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<TestGroupNoInit>(*this);
    }
};

int main() {
    std::cout << "== Wilson GROUP UNIT ==\n";

    auto proxy   = std::make_shared<SpyProxy>();
    auto comp    = std::make_shared<SpyComposer>();
    auto use_mty = std::make_shared<DummyBoolAPI>();
    auto mname   = std::make_shared<DummyStringAPI>();
    auto mpath   = std::make_shared<DummyPathAPI>();
    WilsonGroupAdapterConfig cfg(proxy, comp, use_mty, mname, mpath);

    TestGroupNoInit grp(WGroup::B, cfg, ContributionType::SM);

    auto storage_match = GroupMapper::str(WGroup::B, ScaleType::MATCHING);
    grp.set_matching_storage_block(storage_match);
    assert(grp.get_matching_storage_block() == storage_match);

    auto c7 = std::make_shared<C7>();
    grp.insert({{"C7", c7}});

    grp.init(QCDOrder::LO);

    assert(c7->get_storage_block() == storage_match);
    assert(c7->get_type() == ContributionType::SM);

    assert(!comp->params.empty());

    bool found_LO_target = false;
    for (auto& r : comp->params) {
        ParamId expected_target{storage_match, c7->get_lhaid(QCDOrder::LO)};
        if (r.target == expected_target) {
            ParamId expected_src{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)};
            assert(r.sources.count(expected_src) == 1);
            found_LO_target = true;
        }
    }
    assert(found_LO_target);

    proxy->ret = 7.5;
    auto v = grp.get_matching_coefficient("C7", "LO", ContributionType::SM);
    assert(std::abs(v.real() - 7.5) < 1e-12);

    proxy->ret = 2.25;
    auto vr = grp.get_running_coefficient("C7", "LO", ContributionType::SM, WilsonBasis::B_STANDARD);
    assert(std::abs(vr.real() - 2.25) < 1e-12);

    auto had_block = GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD);
    assert(proxy->last_block == had_block);
    assert(proxy->last_id == c7->id(QCDOrder::LO, ContributionType::SM));

    assert(grp.get_order() == QCDOrder::LO);

    auto gcopy = std::make_shared<TestGroupNoInit>(grp);
    assert(gcopy->count("C7") == 1);
    assert(gcopy->at("C7").get() != grp.at("C7").get());
    assert(gcopy->at("C7")->get_storage_block() == grp.at("C7")->get_storage_block());

    std::cout << "UNIT OK\n";
    return 0;
}
