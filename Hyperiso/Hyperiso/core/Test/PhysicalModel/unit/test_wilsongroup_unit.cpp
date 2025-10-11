#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "WilsonGroup.h"    // CoefficientGroup + WilsonGroupAdapterConfig
#include "Wilson.h"         // WilsonCoefficient
#include "BWilson.h"             // ton coefficient concret
#include "IParamAdapter.h"  // IParameterProxy
#include "IBlockComposer.h" // IBlockComposer
#include "ICoreAPI.h"       // ICoreAPI
#include "Parameter.h"      // ParamId, Parameter, Block
#include "Include.h"        // enums, GroupMapper::str, OrderMapper...

// ---- DUMMIES (ports) -------------------------------------------------

class DummyBoolAPI : public ICoreAPI<bool> {
public: bool get() override { return false; } };

class DummyStringAPI : public ICoreAPI<std::string> {
public: std::string get() override { return "SM"; } };

class DummyPathAPI : public ICoreAPI<fs::path> {
public: fs::path get() override { return fs::path("/dev/null"); } };

// Proxy espion IParameterProxy
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

// Composer espion IBlockComposer
struct ParamRecord {
    ParamId target;
    std::unordered_set<ParamId> sources;
};

class SpyComposer : public IBlockComposer {
public:
    std::vector<ParamRecord> params;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override
    { /* pas utilisé ici */ }

    void compose_parameter(const ParamId& id,
                           const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc&) override
    {
        params.push_back({id, src});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

// Groupe concret minimal pour pouvoir fixer l’id (WGroup)
class TestGroup : public CoefficientGroup {
public:
    using CoefficientGroup::CoefficientGroup;
    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<TestGroup>(*this);
    }
    void set_id(WGroup gid) { this->id = gid; }
    void set_type(ContributionType t) { this->wilson_type = t; }
};

int main() {
    std::cout << "== Wilson GROUP UNIT ==\n";

    auto proxy   = std::make_shared<SpyProxy>();
    auto comp    = std::make_shared<SpyComposer>();
    auto use_mty = std::make_shared<DummyBoolAPI>();
    auto mname   = std::make_shared<DummyStringAPI>();
    auto mpath   = std::make_shared<DummyPathAPI>();

    WilsonGroupAdapterConfig cfg(proxy, comp, use_mty, mname, mpath);

    // Un seul coefficient concret (C7)
    std::map<std::string, std::shared_ptr<WilsonCoefficient>> coeffs;
    auto c7 = std::make_shared<C7>();
    coeffs["C7"] = c7;

    TestGroup grp(coeffs, cfg);
    grp.set_id(WGroup::B);             // le groupe B
    grp.set_type(ContributionType::SM);

    // 1) Claim : storage block bien fixé et type propagé
    auto storage_match = GroupMapper::str(WGroup::B, ScaleType::MATCHING);
    assert(c7->get_storage_block() == storage_match);
    assert(c7->get_type() == ContributionType::SM);

    // 2) init() a appelé compose_parameter au moins pour LO
    assert(!comp->params.empty());
    bool found_LO_target = false;
    for (auto& r : comp->params) {
        if (r.target == ParamId{storage_match, c7->get_lhaid(QCDOrder::LO)}) {
            // sources LO doivent contenir WPARAM_MATCH_SM, (2,1)
            ParamId expected{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)};
            assert(r.sources.count(expected) == 1);
            found_LO_target = true;
        }
    }
    assert(found_LO_target);

    // 3) matching: lecture via proxy
    proxy->ret = 7.5;
    auto v = grp.get_matching_coefficient("C7", "LO", ContributionType::SM);
    assert(std::abs(v.real() - 7.5) < 1e-12);

    // 4) running: lecture via proxy (bloc hadronique)
    proxy->ret = 2.25;
    auto vr = grp.get_running_coefficient("C7", "LO", ContributionType::SM, WilsonBasis::B_STANDARD);
    assert(std::abs(vr.real() - 2.25) < 1e-12);
    auto had_block = GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD);
    assert(proxy->last_block == had_block);
    assert(proxy->last_id == c7->id(QCDOrder::LO, ContributionType::SM));

    // 5) get_order() — par défaut le scan choisit LO si max_order non renseigné
    assert(grp.get_order() == QCDOrder::LO);

    // 6) copie profonde du groupe (clone des WC)
    auto gcopy = std::make_shared<TestGroup>(grp);
    assert(gcopy->count("C7") == 1);
    // pointeurs différents
    assert(gcopy->at("C7").get() != grp.at("C7").get());
    // mêmes propriétés
    assert(gcopy->at("C7")->get_name() == grp.at("C7")->get_name());
    assert(gcopy->at("C7")->get_storage_block() == grp.at("C7")->get_storage_block());
    assert(gcopy->at("C7")->get_type() == grp.at("C7")->get_type());

    std::cout << "✅ UNIT OK\n";
    return 0;
}
