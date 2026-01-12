#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
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

#include "BllDecay.h"
#include "ObsWilsonHelper.h"
#include "ObsPortsConfig.h"
#include "IObsWilsonProxy.h"
#include "Parameter.h"
#include "Math.h"


class FakeCoreBool : public IObsCoreAPI<bool> {
public:
    bool v = false;
    explicit FakeCoreBool(bool b=false) : v(b) {}
    bool get() override { return v; }
};

class SpyObsParameterProxy : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
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

class FakeObsQCDProxy : public IObsQCDProxy {
public:
    double operator()(AlphasConfig) override { return 0.0; }
    double operator()(MassConfig) override { return 0.0; }
    QCDConstants* get_constants() override { return nullptr; }
};

class SpyWilsonFreezer : public IWilsonFreezer<WGroupId> {
public:
    std::vector<WGroupId> froze, unfroze;
    void freeze(WGroupId g) override   { froze.push_back(g); }
    void unfreeze(WGroupId g) override { unfroze.push_back(g); }
};

enum class WKind { M, FM, R, FR };
struct WKey {
    WKind kind; int g; int c; int o; int t;
    bool operator==(const WKey& other) const {
        return kind==other.kind && g==other.g && c==other.c && o==other.o && t==other.t;
    }
};
struct WKeyHash {
    std::size_t operator()(const WKey& k) const noexcept {
        std::size_t h = 1469598103934665603ull;
        auto mix = [&](std::size_t x){ h ^= x + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2); };
        mix((std::size_t)k.kind);
        mix((std::size_t)k.g);
        mix((std::size_t)k.c);
        mix((std::size_t)k.o);
        mix((std::size_t)k.t);
        return h;
    }
};

class SpyObsWilsonProxy : public IObsWilsonProxy {
public:
    bool basis_called = false;
    WilsonBasis last_basis = WilsonBasis::B_STANDARD;
    std::unordered_map<WKey, complex_t, WKeyHash> table;

    void set_value(WKind kind, WGroup g, WCoef c, QCDOrder o, ContributionType t, complex_t v) {
        table[WKey{kind, (int)g, (int)c, (int)o, (int)t}] = v;
    }
    complex_t lookup(WKind kind, WGroup g, WCoef c, QCDOrder o, ContributionType t) {
        auto it = table.find(WKey{kind, (int)g, (int)c, (int)o, (int)t});
        if (it == table.end()) return complex_t{0.0, 0.0};
        return it->second;
    }

    complex_t getM (WGroup g, WCoef c, QCDOrder o, ContributionType t) override { return lookup(WKind::M,  g, c, o, t); }
    complex_t getFM(WGroup g, WCoef c, QCDOrder o, ContributionType t) override { return lookup(WKind::FM, g, c, o, t); }
    complex_t getR (WGroup g, WCoef c, QCDOrder o, ContributionType t) override { return lookup(WKind::R,  g, c, o, t); }
    complex_t getFR(WGroup g, WCoef c, QCDOrder o, ContributionType t) override { return lookup(WKind::FR, g, c, o, t); }

    std::map<QCDOrder, complex_t> getSM(WGroup, WCoef, ContributionType) override { return {}; }
    std::map<QCDOrder, complex_t> getSR(WGroup, WCoef, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAM (WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAR (WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAFM(WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAFR(WGroup, QCDOrder, ContributionType) override { return {}; }

    std::shared_ptr<IObsWilsonBuilder> get_builder() override { return nullptr; }
    std::unordered_set<WilsonBasis> get_bases(WGroupId) override { return {WilsonBasis::B_STANDARD}; }

    void set_basis(WilsonBasis b) override { basis_called = true; last_basis = b; }
};

class SpyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    int build_calls = 0;
    std::shared_ptr<SpyObsWilsonProxy> proxy = std::make_shared<SpyObsWilsonProxy>();

    void build(std::shared_ptr<AbstractConfig>) override { build_calls++; }
    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return proxy; }
};

static inline scalar_t R(double x) { return scalar_t{x, 0.0}; }
static inline scalar_t C(double re, double im) { return scalar_t{re, im}; }

static double expected_BR_Bs_mumu(const SpyObsParameterProxy& pp,
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

    const double qcd_52 = pp.values.at(ParamId{ParameterType::SM, "QCD", {5,2}}).real();
    const double m_s    = pp.values.at(ParamId{ParameterType::SM, "MASS", 3}).real();
    const double r_s    = m_Bs / (qcd_52 + m_s);

    const complex_t C10_SM = wp.getFR(WGroup::B,       WCoef::C10,  order, ContributionType::SM);
    const complex_t C10    = wp.getFR(WGroup::B,       WCoef::C10,  order, ContributionType::TOTAL);
    const complex_t CQ1    = wp.getFR(WGroup::BScalar, WCoef::CQ1,  order, ContributionType::TOTAL);
    const complex_t CQ2    = wp.getFR(WGroup::BScalar, WCoef::CQ2,  order, ContributionType::TOTAL);
    const complex_t CP10   = wp.getFR(WGroup::BPrime,  WCoef::CP10, order, ContributionType::TOTAL);
    const complex_t CPQ1   = wp.getFR(WGroup::BPrime,  WCoef::CPQ1, order, ContributionType::TOTAL);
    const complex_t CPQ2   = wp.getFR(WGroup::BPrime,  WCoef::CPQ2, order, ContributionType::TOTAL);

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
    std::cout << "== BllDecay INTEGRATION ==\n";

    ObsWilsonHelper(true);

    auto builder_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> builder = builder_spy;

    auto p_spy = std::make_shared<SpyObsParameterProxy>();
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_sm   = p_spy;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_flav = p_spy;

    auto qcd = std::make_shared<FakeObsQCDProxy>();
    auto use_marty = std::make_shared<FakeCoreBool>(false);

    auto freezer_spy = std::make_shared<SpyWilsonFreezer>();
    std::shared_ptr<IWilsonFreezer<WGroupId>> freezer = freezer_spy;

    ObservablePortsConfig ports(builder, p_sm, p_flav, qcd, use_marty, freezer);

    p_spy->values[ParamId{ParameterType::SM,     "SMINPUTS", 2}]       = R(1.0e-5);
    p_spy->values[ParamId{ParameterType::SM,     "EW", {1,2}}]         = R(1.0/137.0);
    p_spy->values[ParamId{ParameterType::SM,     "MASS", 13}]          = R(0.105);
    p_spy->values[ParamId{ParameterType::FLAVOR, "FMASS", 531}]        = R(5.367);
    p_spy->values[ParamId{ParameterType::FLAVOR, "FCONST", {531,1}}]   = R(0.230);
    p_spy->values[ParamId{ParameterType::FLAVOR, "FLIFE", 531}]        = R(HBAR * 1.6);

    p_spy->values[ParamId{ParameterType::SM, "VCKM", {2,2}}]           = C(1.0, 0.0);
    p_spy->values[ParamId{ParameterType::SM, "VCKM", {2,1}}]           = C(1.0, 0.0);

    p_spy->values[ParamId{ParameterType::DECAY, "B_ll", 2}]            = R(1.0);

    p_spy->values[ParamId{ParameterType::SM, "QCD", {5,2}}]            = R(4.2);
    p_spy->values[ParamId{ParameterType::SM, "MASS", 3}]               = R(0.09);

    auto& W = *builder_spy->proxy;
    const QCDOrder ord = QCDOrder::LO;

    W.set_value(WKind::FR, WGroup::B,       WCoef::C10,  ord, ContributionType::SM,    complex_t{-4.0, 0.0});
    W.set_value(WKind::FR, WGroup::B,       WCoef::C10,  ord, ContributionType::TOTAL, complex_t{-3.5, 0.0});
    W.set_value(WKind::FR, WGroup::BScalar, WCoef::CQ1,  ord, ContributionType::TOTAL, complex_t{ 0.2, 0.0});
    W.set_value(WKind::FR, WGroup::BScalar, WCoef::CQ2,  ord, ContributionType::TOTAL, complex_t{-0.1, 0.0});
    W.set_value(WKind::FR, WGroup::BPrime,  WCoef::CP10, ord, ContributionType::TOTAL, complex_t{ 0.1, 0.0});
    W.set_value(WKind::FR, WGroup::BPrime,  WCoef::CPQ1, ord, ContributionType::TOTAL, complex_t{ 0.05,0.0});
    W.set_value(WKind::FR, WGroup::BPrime,  WCoef::CPQ2, ord, ContributionType::TOTAL, complex_t{-0.02,0.0});

    BllDecay decay(QCDOrder::NONE, 160.0, 4.8, ports);
    decay.set_order(QCDOrder::LO);

    decay.enable();
    assert(builder_spy->build_calls == 1);
    assert(W.basis_called);
    assert(W.last_basis == WilsonBasis::B_STANDARD);

    {
        auto out = decay.compute_observable(Observables::BR_BS_MUMU);
        assert(out.size() == 1);

        const double expected = expected_BR_Bs_mumu(*p_spy, W, QCDOrder::LO);
        assert(std::abs(out[0].value - expected) / (std::abs(expected) + 1e-30) < 1e-12);
    }

    {
        const double out1 = decay.compute_observable(Observables::BR_BS_MUMU)[0].value;

        p_spy->values[ParamId{ParameterType::FLAVOR, "FCONST", {531,1}}] = R(0.240);

        decay.disable();
        decay.enable();

        assert(builder_spy->build_calls == 1);

        const double out2 = decay.compute_observable(Observables::BR_BS_MUMU)[0].value;
        assert(out1 != out2);
    }

    {
        ObsWilsonHelper(true);
        use_marty->v = true;

        BllDecay decay2(QCDOrder::NONE, 160.0, 4.8, ports);
        decay2.set_order(QCDOrder::NLO);
        decay2.enable();

        auto out = decay2.compute_observable(Observables::BR_BS_MUMU);
        assert(!out.empty());
    }

    std::cout << "INTEGRATION OK\n";
    return 0;
}
