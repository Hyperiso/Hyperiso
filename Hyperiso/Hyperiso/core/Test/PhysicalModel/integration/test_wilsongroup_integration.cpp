#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <map>

#include "WilsonGroup.h"
#include "Wilson.h"
#include "BWilsonGroup.h"
#include "IParamAdapter.h"
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "Parameter.h"
#include "Include.h"

// ---- DUMMIES ports ----------------------------------------------------
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
    struct R { ParamId target; std::unordered_set<ParamId> sources; QCDOrder order; };
    std::vector<R> calls;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& id,
                           const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc&) override {
        // On ne connaît pas l’ordre ici, mais on l’infère à partir de LhaID (3e champ: order)
        QCDOrder ord = QCDOrder::LO;
        // if (id.id.vals.size() >= 3) {
        //     // LhaID(cte) -> dépend de ta struct; fallback LO si indispo
        //     // Ici on s’en passe : on stocke LO par défaut
        // }
        calls.push_back({id, src, ord});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

// ---- Coeff factice avec max_order=NLO ---------------------------------
class DummyCoeff : public WilsonCoefficient {
public:
    DummyCoeff() : WilsonCoefficient("DUM", "X_MATCH") {
        // on indique NLO comme max_order
        this->max_order = QCDOrder::NLO;

        // LO: 1 source
        matching_info[QCDOrder::LO] = {
            { {ParameterType::WILSON, "DUMSRC", LhaID(1)} },
            [](const auto&){ return 1.0; },
            LhaID(111,222,0,0)
        };
        // NLO: 2 sources
        matching_info[QCDOrder::NLO] = {
            {
                {ParameterType::WILSON, "DUMSRC", LhaID(2)},
                {ParameterType::SM,     "MASS",   LhaID(24)}
            },
            [](const auto&){ return 2.0; },
            LhaID(111,222,1,0)
        };
    }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<DummyCoeff>(*this);
    }
};

// Groupe concret
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
    std::cout << "== Wilson GROUP INTEGRATION ==\n";

    auto proxy   = std::make_shared<SpyProxy>();
    auto comp    = std::make_shared<SpyComposer>();
    auto use_mty = std::make_shared<DummyBoolAPI>();
    auto mname   = std::make_shared<DummyStringAPI>();
    auto mpath   = std::make_shared<DummyPathAPI>();

    WilsonGroupAdapterConfig cfg(proxy, comp, use_mty, mname, mpath);

    std::map<std::string, std::shared_ptr<WilsonCoefficient>> coeffs;
    coeffs["C7"]  = std::make_shared<C7>();       // max_order par défaut -> LO
    coeffs["DUM"] = std::make_shared<DummyCoeff>(); // max_order NLO

    TestGroup grp(coeffs, cfg);
    grp.set_id(WGroup::B);
    grp.set_type(ContributionType::SM);

    // 1) Composer appelé au moins 1 (C7 LO) + 2 (DUM LO+NLO) => >= 3 appels
    assert(comp->calls.size() >= 3);

    // 2) Vérifie qu’on a bien la cible de C7 LO dans le bloc MATCHING du groupe
    const auto match_blk = GroupMapper::str(WGroup::B, ScaleType::MATCHING);
    bool saw_c7_lo = false;
    for (const auto& r : comp->calls) {
        if (r.target == ParamId{match_blk, coeffs["C7"]->get_lhaid(QCDOrder::LO)}) {
            // doit contenir WPARAM_MATCH_SM, (2,1)
            ParamId need{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)};
            assert(r.sources.count(need) == 1);
            saw_c7_lo = true;
        }
    }
    assert(saw_c7_lo);

    // 3) Vérifie qu’on a bien les deux ordres pour DUM
    bool saw_dum_lo  = false;
    bool saw_dum_nlo = false;
    for (const auto& r : comp->calls) {
        if (r.target == ParamId{match_blk, coeffs["DUM"]->get_lhaid(QCDOrder::LO)})  saw_dum_lo  = true;
        if (r.target == ParamId{match_blk, coeffs["DUM"]->get_lhaid(QCDOrder::NLO)}) saw_dum_nlo = true;
    }
    assert(saw_dum_lo && saw_dum_nlo);

    // 4) get_running_coefficient : vérifie le bloc hadronique + id
    proxy->ret = 9.0;
    auto vr = grp.get_running_coefficient("DUM", "NLO", ContributionType::SM, WilsonBasis::B_STANDARD);
    assert(std::abs(vr.real() - 9.0) < 1e-12);
    auto had_blk = GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD);
    assert(proxy->last_block == had_blk);
    assert(proxy->last_id == coeffs["DUM"]->id(QCDOrder::NLO, ContributionType::SM));

    // 5) L’ordre du groupe vaut le max des coeffs (ici NLO grâce à DUM)
    assert(grp.get_order() == QCDOrder::NLO);

    std::cout << "✅ INTEGRATION OK\n";
    return 0;
}
