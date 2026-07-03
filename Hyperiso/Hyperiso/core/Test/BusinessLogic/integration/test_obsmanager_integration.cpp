#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <cmath>

#include "Include.h"

#ifdef LOG_DEBUG
#undef LOG_DEBUG
#endif
#ifdef LOG_INFO
#undef LOG_INFO
#endif
#ifdef LOG_WARN
#undef LOG_WARN
#endif
#ifdef LOG_ERROR
#undef LOG_ERROR
#endif
#define LOG_DEBUG(...) do {} while(0)
#define LOG_INFO(...)  do {} while(0)
#define LOG_WARN(...)  do {} while(0)
#define LOG_ERROR(...) do {} while(0)

#include "ObsManager.h"
#include "ObsWilsonHelper.h"
#include "IObsWilsonProxy.h"
#include "IObsWilsonBuilder.h"
#include "IObsParameterProxy.h"
#include "IObsCoreAPI.h"
#include "IObsQCDProxy.h"
#include "IWilsonFreezer.h"
#include "Math.h"


class FakeCoreBool : public IObsCoreAPI<bool> {
public:
    bool v=false;
    explicit FakeCoreBool(bool b=false): v(b) {}
    bool get() override { return v; }
};

class FakeObsQCDProxy : public IObsQCDProxy {
public:
    double operator()(AlphasConfig) override { return 0.0; }
    double operator()(MassConfig) override { return 0.0; }
    QCDConstants* get_constants() override { return nullptr; }
};

class SpyFreezer : public IWilsonFreezer<WGroupId> {
public:
    void freeze(WGroupId) override {}
    void unfreeze(WGroupId) override {}
};

class SpyParamProxy : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    std::unordered_map<ParamId, scalar_t> values;

    scalar_t operator()(const ParamId& pid, DataType) override {
        auto it = values.find(pid);
        if (it == values.end()) return scalar_t{0.0, 0.0};
        return it->second;
    }

    scalar_t operator()(const std::string&, const LhaID&, DataType) const override {
        return scalar_t{0.0, 0.0};
    }

    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const override {
        return std::make_shared<Parameter>(pid, 0.0, 0.0, 0.0);
    }
};

struct WKey {
    int g, c, o, t;
    bool operator==(const WKey& other) const { return g==other.g && c==other.c && o==other.o && t==other.t; }
};
struct WKeyHash {
    std::size_t operator()(const WKey& k) const noexcept {
        std::size_t h=0;
        h ^= (std::size_t)k.g + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.c + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.o + 0x9e37 + (h<<6) + (h>>2);
        h ^= (std::size_t)k.t + 0x9e37 + (h<<6) + (h>>2);
        return h;
    }
};

class SpyObsWilsonProxy : public IObsWilsonProxy {
public:
    bool basis_called=false;
    WilsonBasis last_basis=WilsonBasis::B_STANDARD;

    std::unordered_map<WKey, complex_t, WKeyHash> fr;

    complex_t getM (WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getFM(WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getR (WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }

    complex_t getFR(WGroup g, WCoef c, QCDOrder o, ContributionType t) override {
        auto it = fr.find(WKey{(int)g,(int)c,(int)o,(int)t});
        if (it == fr.end()) return {0,0};
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

    void set_basis(WilsonBasis b) override { basis_called=true; last_basis=b; }
};

class SpyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    int build_calls=0;
    std::shared_ptr<SpyObsWilsonProxy> proxy = std::make_shared<SpyObsWilsonProxy>();

    void build(std::shared_ptr<AbstractConfig>) override { build_calls++; }
    void add_custom_group(const CustomWilsonGroupConfig&) override {}
    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return proxy; }
};

static inline scalar_t R(double x){ return scalar_t{x,0.0}; }
static inline scalar_t C(double re,double im){ return scalar_t{re,im}; }

static double expected_BR_Bs_mumu(const SpyParamProxy& pp,
                                 SpyObsWilsonProxy& wp,
                                 QCDOrder order)
{
    const double G_F      = pp.values.at(ParamId{ParameterType::SM, "SMINPUTS", 2}).real();
    const double alpha_em = pp.values.at(ParamId{ParameterType::SM, "EW", {1,2}}).real();
    const double m_mu     = pp.values.at(ParamId{ParameterType::SM, "MASS", 13}).real();
    const double m_Bs     = pp.values.at(ParamId{ParameterType::FLAVOR, "FMASS", 531}).real();
    const double f_Bs     = pp.values.at(ParamId{ParameterType::FLAVOR, "FCONST", {531,1}}).real();
    const double tau_Bs   = pp.values.at(ParamId{ParameterType::FLAVOR, "FLIFE", 531}).real() / HBAR;

    const complex_t V22 = pp.values.at(ParamId{ParameterType::SM, "VCKM", {2,2}});
    const complex_t V21 = pp.values.at(ParamId{ParameterType::SM, "VCKM", {2,1}});
    const complex_t lambda_s = V22 * std::conj(V21);

    const double eta_BBS = pp.values.at(ParamId{ParameterType::DECAY, "B_ll", 2}).real();

    const double x_s = m_mu / m_Bs;
    const double beta_s = std::sqrt(1.0 - 4.0 * x_s * x_s);

    const double mb = pp.values.at(ParamId{ParameterType::SM, "QCD", {5,2}}).real();
    const double ms = pp.values.at(ParamId{ParameterType::SM, "MASS", 3}).real();
    const double r_s = m_Bs / (mb + ms);

    const complex_t C10_SM = wp.getFR(WGroup::B,      WCoef::C10,  order, ContributionType::SM);
    const complex_t C10    = wp.getFR(WGroup::B,      WCoef::C10,  order, ContributionType::TOTAL);
    const complex_t CQ1    = wp.getFR(WGroup::BScalar,WCoef::CQ1_MU,  order, ContributionType::TOTAL);
    const complex_t CQ2    = wp.getFR(WGroup::BScalar,WCoef::CQ2_MU,  order, ContributionType::TOTAL);

    const complex_t CP10   = wp.getFR(WGroup::BPrime, WCoef::CP10, order, ContributionType::TOTAL);
    const complex_t CPQ1   = wp.getFR(WGroup::BPrime, WCoef::CPQ1_MU, order, ContributionType::TOTAL);
    const complex_t CPQ2   = wp.getFR(WGroup::BPrime, WCoef::CPQ2_MU, order, ContributionType::TOTAL);

    const complex_t C10_m = C10 - CP10;
    const complex_t CQ1_m = CQ1 - CPQ1;
    const complex_t CQ2_m = CQ2 - CPQ2;

    const double pref = std::pow(G_F * alpha_em, 2) / (64.0 * std::pow(M_PI, 3)) * eta_BBS;
    const double term1 = std::pow(beta_s * std::abs(r_s * CQ1_m), 2);
    const double term2 = std::pow(std::abs(r_s * CQ2_m + 2.0 * x_s * C10_m), 2);

    return pref * std::pow(f_Bs * std::abs(lambda_s), 2)
                * std::pow(m_Bs, 3)
                * tau_Bs
                * beta_s
                * (term1 + term2);
}

int main() {
    std::cout << "== ObsManager INTEGRATION ==\n";

    ObsWilsonHelper(true);

    auto wb_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> wb = wb_spy;

    auto pp = std::make_shared<SpyParamProxy>();
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_sm   = pp;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_flav = pp;

    auto qcd = std::make_shared<FakeObsQCDProxy>();
    auto use_marty = std::make_shared<FakeCoreBool>(false);
    auto freezer = std::make_shared<SpyFreezer>();

    pp->values[ParamId{ParameterType::SM, "MASS", 24}] = R(80.0);
    pp->values[ParamId{ParameterType::SM, "QCD", LhaID(5,2)}] = R(4.2);

    pp->values[ParamId{ParameterType::SM,     "SMINPUTS", 2}]      = R(1.0e-5);
    pp->values[ParamId{ParameterType::SM,     "EW", {1,2}}]        = R(1.0/137.0);
    pp->values[ParamId{ParameterType::SM,     "MASS", 13}]         = R(0.105);

    pp->values[ParamId{ParameterType::FLAVOR, "FMASS", 511}]       = R(5.279);
    pp->values[ParamId{ParameterType::FLAVOR, "FMASS", 531}]       = R(5.367);

    pp->values[ParamId{ParameterType::FLAVOR, "FCONST", {511,1}}]  = R(0.190);
    pp->values[ParamId{ParameterType::FLAVOR, "FCONST", {531,1}}]  = R(0.230);

    pp->values[ParamId{ParameterType::FLAVOR, "FLIFE", 511}]       = R(HBAR * 1.5);
    pp->values[ParamId{ParameterType::FLAVOR, "FLIFE", 531}]       = R(HBAR * 1.6);

    pp->values[ParamId{ParameterType::SM, "VCKM", {2,2}}]          = C(1.0, 0.0);
    pp->values[ParamId{ParameterType::SM, "VCKM", {2,0}}]          = C(1.0, 0.0);
    pp->values[ParamId{ParameterType::SM, "VCKM", {2,1}}]          = C(1.0, 0.0);

    pp->values[ParamId{ParameterType::DECAY, "B_ll", 1}]           = R(0.0); // ys
    pp->values[ParamId{ParameterType::DECAY, "B_ll", 2}]           = R(1.0); // eta_BBS

    pp->values[ParamId{ParameterType::SM, "MASS", 1}]              = R(0.005);
    pp->values[ParamId{ParameterType::SM, "MASS", 3}]              = R(0.09);
    pp->values[ParamId{ParameterType::SM, "QCD", {5,2}}]           = R(4.2);

    ObservablePortsConfig ports(wb, p_sm, p_flav, qcd, use_marty, freezer);

    auto& W = *wb_spy->proxy;
    const QCDOrder ord = QCDOrder::LO;

    W.fr[WKey{(int)WGroup::B,      (int)WCoef::C10,  (int)ord, (int)ContributionType::SM}]    = complex_t{-4.0, 0.0};
    W.fr[WKey{(int)WGroup::B,      (int)WCoef::C10,  (int)ord, (int)ContributionType::TOTAL}] = complex_t{-3.5, 0.0};
    W.fr[WKey{(int)WGroup::BScalar,(int)WCoef::CQ1_MU,  (int)ord, (int)ContributionType::TOTAL}] = complex_t{ 0.2, 0.0};
    W.fr[WKey{(int)WGroup::BScalar,(int)WCoef::CQ2_MU,  (int)ord, (int)ContributionType::TOTAL}] = complex_t{-0.1, 0.0};
    W.fr[WKey{(int)WGroup::BPrime, (int)WCoef::CP10, (int)ord, (int)ContributionType::TOTAL}] = complex_t{ 0.1, 0.0};
    W.fr[WKey{(int)WGroup::BPrime, (int)WCoef::CPQ1_MU, (int)ord, (int)ContributionType::TOTAL}] = complex_t{ 0.05,0.0};
    W.fr[WKey{(int)WGroup::BPrime, (int)WCoef::CPQ2_MU, (int)ord, (int)ContributionType::TOTAL}] = complex_t{-0.02,0.0};

    ObsManager mgr(ports, false);
    auto did = DecayMapper::to_id(Decays::B__l_l);
    mgr.add_custom_decay(did, std::make_shared<BllDecay>(QCDOrder::NONE, 160.0, 4.8, ports));

    mgr.add_obs(Observables::BR_BS_MUMU, QCDOrder::LO, false);

    auto out = mgr.evaluate(Observables::BR_BS_MUMU);
    assert(out.size() == 1);

    const double expected = expected_BR_Bs_mumu(*pp, W, QCDOrder::LO);
    assert(std::abs(out[0].value - expected) / (std::abs(expected) + 1e-30) < 1e-12);

    {
        auto all = mgr.evaluate_all();
        const auto oid = ObservableMapper::to_id(Observables::BR_BS_MUMU);
        assert(all.count(oid) == 1);
        assert(all[oid].size() == 1);
    }

    std::cout << "INTEGRATION OK\n";
    return 0;
}
