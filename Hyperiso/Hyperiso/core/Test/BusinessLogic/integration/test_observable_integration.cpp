#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <map>
#include <cmath>

#include "Include.h"

#ifdef LOG_DEBUG
#undef LOG_DEBUG
#endif
#ifdef LOG_WARN
#undef LOG_WARN
#endif
#ifdef LOG_ERROR
#undef LOG_ERROR
#endif
#define LOG_DEBUG(...) do {} while(0)
#define LOG_WARN(...)  do {} while(0)
#define LOG_ERROR(...) do {} while(0)

#include "Observable.h"
#include "BllDecay.h"
#include "ObsWilsonHelper.h"
#include "ObsPortsConfig.h"
#include "IObsWilsonProxy.h"
#include "Parameter.h"
#include "Math.h"


class SpyObsParamProxy : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    std::map<std::pair<ParamId, DataType>, scalar_t> table;

    scalar_t operator()(const ParamId& pid, DataType dt) override {
        auto it = table.find({pid, dt});
        if (it == table.end()) return scalar_t{0.0, 0.0};
        return it->second;
    }
    scalar_t operator()(const std::string&, const LhaID&, DataType) const override { return scalar_t{0.0, 0.0}; }
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const override {
        return std::make_shared<Parameter>(pid, 0.0, 0.0, 0.0);
    }
};

class FakeCoreBool : public IObsCoreAPI<bool> {
public:
    bool v=false;
    explicit FakeCoreBool(bool b=false) : v(b) {}
    bool get() override { return v; }
};

class SpyObsParameterProxySM : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    std::unordered_map<ParamId, scalar_t> values;

    scalar_t operator()(const ParamId& pid, DataType) override {
        auto it = values.find(pid);
        if (it == values.end()) return scalar_t{0.0, 0.0};
        return it->second;
    }
    scalar_t operator()(const std::string&, const LhaID&, DataType) const override { return scalar_t{0.0, 0.0}; }
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const override {
        return std::make_shared<Parameter>(pid, 0.0, 0.0, 0.0);
    }
};

class FakeObsQCDProxy : public IObsQCDProxy {
public:
    double operator()(AlphasConfig) override { return 0.0; }
    double operator()(MassConfig) override { return 0.0; }
    QCDConstants* get_constants() override { return nullptr; }
};

class SpyWilsonFreezer : public IWilsonFreezer<WGroupId> {
public:
    void freeze(WGroupId) override {}
    void unfreeze(WGroupId) override {}
};

enum class WKind { M, FM, R, FR };
struct WKey { WKind kind; int g; int c; int o; int t;
    bool operator==(const WKey& other) const {
        return kind==other.kind && g==other.g && c==other.c && o==other.o && t==other.t;
    }
};
struct WKeyHash {
    std::size_t operator()(const WKey& k) const noexcept {
        std::size_t h=0;
        h ^= (std::size_t)k.kind + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.g    + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.c    + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.o    + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.t    + 0x9e37 + (h<<6) + (h>>2);
        return h;
    }
};

class SpyObsWilsonProxy : public IObsWilsonProxy {
public:
    std::unordered_map<WKey, complex_t, WKeyHash> table;
    void set_basis(WilsonBasis) override {}

    complex_t getM (WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getFM(WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getR (WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getFR(WGroup g, WCoef c, QCDOrder o, ContributionType t) override {
        auto it = table.find(WKey{WKind::FR,(int)g,(int)c,(int)o,(int)t});
        if (it==table.end()) return {0,0};
        return it->second;
    }


    static WGroup enum_group_or_default(WGroupId gid) {
        auto g = GroupMapper::enum_of(gid);
        return g.value_or(WGroup::B);
    }

    static WCoef enum_coef_or_default(WCoefId cid) {
        auto c = WCoefMapper::enum_of(cid);
        return c.value_or(WCoef::C1);
    }

    complex_t getM (WGroupId g, WCoefId c, QCDOrder o, ContributionType t) override { return getM (enum_group_or_default(g), enum_coef_or_default(c), o, t); }
    complex_t getFM(WGroupId g, WCoefId c, QCDOrder o, ContributionType t) override { return getFM(enum_group_or_default(g), enum_coef_or_default(c), o, t); }
    complex_t getR (WGroupId g, WCoefId c, QCDOrder o, ContributionType t) override { return getR (enum_group_or_default(g), enum_coef_or_default(c), o, t); }
    complex_t getFR(WGroupId g, WCoefId c, QCDOrder o, ContributionType t) override { return getFR(enum_group_or_default(g), enum_coef_or_default(c), o, t); }

    std::map<QCDOrder, complex_t> getSM(WGroupId g, WCoefId c, ContributionType t) override { return getSM(enum_group_or_default(g), enum_coef_or_default(c), t); }
    std::map<QCDOrder, complex_t> getSR(WGroupId g, WCoefId c, ContributionType t) override { return getSR(enum_group_or_default(g), enum_coef_or_default(c), t); }

    std::map<QCDOrder, complex_t> getSM(WGroup, WCoef, ContributionType) override { return {}; }
    std::map<QCDOrder, complex_t> getSR(WGroup, WCoef, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAM (WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAR (WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAFM(WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAFR(WGroup, QCDOrder, ContributionType) override { return {}; }

    std::shared_ptr<IObsWilsonBuilder> get_builder() override { return nullptr; }
    std::unordered_set<WilsonBasis> get_bases(WGroupId) override { return {WilsonBasis::B_STANDARD}; }
};

class SpyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    std::shared_ptr<SpyObsWilsonProxy> proxy = std::make_shared<SpyObsWilsonProxy>();
    void build(std::shared_ptr<AbstractConfig>) override {}
    void add_custom_group(const CustomWilsonGroupConfig&) override {}
    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return proxy; }
};

static inline scalar_t R(double x){ return scalar_t{x,0.0}; }
static inline scalar_t C(double re,double im){ return scalar_t{re,im}; }

int main() {
    std::cout << "== Observable INTEGRATION ==\n";

    ObsWilsonHelper(true);

    auto wb_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> wb = wb_spy;

    auto sm_spy = std::make_shared<SpyObsParameterProxySM>();
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_sm = sm_spy;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_flav = sm_spy;

    auto qcd = std::make_shared<FakeObsQCDProxy>();
    auto use_marty = std::make_shared<FakeCoreBool>(false);
    auto freezer = std::make_shared<SpyWilsonFreezer>();

    ObservablePortsConfig ports(wb, p_sm, p_flav, qcd, use_marty, freezer);

    sm_spy->values[ParamId{ParameterType::SM,     "SMINPUTS", 2}]      = R(1.0e-5);
    sm_spy->values[ParamId{ParameterType::SM,     "EW", {1,2}}]        = R(1.0/137.0);
    sm_spy->values[ParamId{ParameterType::SM,     "MASS", 13}]         = R(0.105);
    sm_spy->values[ParamId{ParameterType::FLAVOR, "FMASS", 531}]       = R(5.367);
    sm_spy->values[ParamId{ParameterType::FLAVOR, "FCONST", {531,1}}]  = R(0.230);
    sm_spy->values[ParamId{ParameterType::FLAVOR, "FLIFE", 531}]       = R(HBAR * 1.6);
    sm_spy->values[ParamId{ParameterType::SM, "VCKM", {2,2}}]          = C(1.0, 0.0);
    sm_spy->values[ParamId{ParameterType::SM, "VCKM", {2,1}}]          = C(1.0, 0.0);
    sm_spy->values[ParamId{ParameterType::DECAY, "B_ll", 1}]           = R(0.0);
    sm_spy->values[ParamId{ParameterType::DECAY, "B_ll", 2}]           = R(1.0);
    sm_spy->values[ParamId{ParameterType::SM, "QCD", {5,2}}]           = R(4.2);
    sm_spy->values[ParamId{ParameterType::SM, "MASS", 3}]              = R(0.09);

    auto& W = *wb_spy->proxy;
    const QCDOrder ord = QCDOrder::LO;
    W.table[WKey{WKind::FR,(int)WGroup::B,(int)WCoef::C10,(int)ord,(int)ContributionType::SM}]    = complex_t{-4.0,0.0};
    W.table[WKey{WKind::FR,(int)WGroup::B,(int)WCoef::C10,(int)ord,(int)ContributionType::TOTAL}] = complex_t{-3.5,0.0};
    W.table[WKey{WKind::FR,(int)WGroup::BScalar,(int)WCoef::CQ1_MU,(int)ord,(int)ContributionType::TOTAL}] = complex_t{0.2,0.0};
    W.table[WKey{WKind::FR,(int)WGroup::BScalar,(int)WCoef::CQ2_MU,(int)ord,(int)ContributionType::TOTAL}] = complex_t{-0.1,0.0};
    W.table[WKey{WKind::FR,(int)WGroup::BPrime,(int)WCoef::CP10,(int)ord,(int)ContributionType::TOTAL}] = complex_t{0.1,0.0};
    W.table[WKey{WKind::FR,(int)WGroup::BPrime,(int)WCoef::CPQ1_MU,(int)ord,(int)ContributionType::TOTAL}] = complex_t{0.05,0.0};
    W.table[WKey{WKind::FR,(int)WGroup::BPrime,(int)WCoef::CPQ2_MU,(int)ord,(int)ContributionType::TOTAL}] = complex_t{-0.02,0.0};

    auto decay = std::make_shared<BllDecay>(QCDOrder::NONE, 160.0, 4.8, ports);
    decay->set_order(QCDOrder::LO);

    auto obspp = std::make_shared<SpyObsParamProxy>();
    const ObservableId oid = ObservableMapper::to_id(Observables::BR_BS_MUMU);
    const auto flha = ObservableMapper::binned_flha_of({oid, {0.,0.}}).value();
    ParamId pid_obs(ParameterType::OBSERVABLE, "FOBS_DEFAULT", flha);

    obspp->table[{pid_obs, DataType::VALUE}] = scalar_t{2.34, 0.0};
    obspp->table[{pid_obs, UncertaintyTypeMapper::d_type(UncertaintyType::COMBINED)}] = scalar_t{0.10,0.0};

    Observable obs(oid, decay, obspp);

    {
        auto v = obs.get_exp_val();
        assert(std::abs(v.real() - 2.34) < 1e-12);
        auto u = obs.get_exp_uncertainty();
        assert(std::abs(u.real() - 0.10) < 1e-12);
    }

    {
        // Observable::compute() does not enable the decay by itself.
        // The direct Observable path therefore enables the decay explicitly first.
        decay->enable();

        auto out = obs.compute();
        assert(out.size() == 1);
        assert(out[0].id == oid);
        assert(std::isfinite(out[0].value));
    }

    std::cout << "INTEGRATION OK\n";
    return 0;
}
