#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <map>
#include <any>

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

#include "DecayParent.h"
#include "ObsWilsonHelper.h"
#include "ObsPortsConfig.h"
#include "IObsCoreAPI.h"
#include "IObsParameterProxy.h"
#include "IObsQCDProxy.h"
#include "IWilsonFreezer.h"
#include "IObsWilsonProxy.h"
#include "Configs.h"
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
    mutable std::vector<ParamId> calls_pid;

    scalar_t operator()(const ParamId& pid, DataType) override {
        calls_pid.push_back(pid);
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
    std::vector<WGroupId> froze;
    std::vector<WGroupId> unfroze;

    void freeze(WGroupId g) override   { froze.push_back(g); }
    void unfreeze(WGroupId g) override { unfroze.push_back(g); }

    void clear() { froze.clear(); unfroze.clear(); }
};


enum class WKind { M, FM, R, FR };

struct WKey {
    WKind kind;
    int g;
    int c;
    int o;
    int t;
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

    std::unordered_map<WGroupId, std::unordered_set<WilsonBasis>> bases;

    void set_value(WKind kind, WGroup g, WCoef c, QCDOrder o, ContributionType t, complex_t v) {
        table[WKey{kind, (int)g, (int)c, (int)o, (int)t}] = v;
    }

    complex_t lookup(WKind kind, WGroup g, WCoef c, QCDOrder o, ContributionType t) {
        auto it = table.find(WKey{kind, (int)g, (int)c, (int)o, (int)t});
        if (it == table.end()) return complex_t{0.0, 0.0};
        return it->second;
    }

    // --- IObsWilsonProxy ---
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

    std::unordered_set<WilsonBasis> get_bases(WGroupId gid) override {
        auto it = bases.find(gid);
        if (it == bases.end()) return {WilsonBasis::B_STANDARD};
        return it->second;
    }

    void set_basis(WilsonBasis b) override {
        basis_called = true;
        last_basis = b;
    }
};

class SpyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    int build_calls = 0;
    std::unordered_set<WGroupId> last_groups;
    QCDOrder last_order = QCDOrder::NONE;

    std::shared_ptr<SpyObsWilsonProxy> proxy = std::make_shared<SpyObsWilsonProxy>();

    void build(std::shared_ptr<AbstractConfig> cfg) override {
        build_calls++;
        last_groups.clear();

        auto wc = std::dynamic_pointer_cast<WilsonBuildConfig>(cfg);
        if (wc) {
            last_groups = wc->groups;
            last_order  = wc->order;
        }
    }

    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return proxy; }
};

class DummyDecay : public DecayParent {
public:
    bool load_called = false;

    DummyDecay(DecayId id,
               double matching_scale,
               double hadronic_scale,
               QCDOrder order,
               ObservablePortsConfig& ports)
        : DecayParent(id, matching_scale, hadronic_scale, order, ports)
    {
        w_config.groups = { GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime) };
        max_order = QCDOrder::NLO;
    }

    void load_params() override {
        load_called = true;
        (void)(*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    }

    std::vector<ObservableValue> compute_observable(Observables) override { return {}; }
    std::vector<ObservableValue> compute_observable(ObservableId) override { return {}; }
    void set_config(std::any) override {}

    QCDOrder current_order() const { return w_config.order; }
    bool is_enabled() const { return enabled; }
};

static std::unordered_set<WGroupId> set_from(std::initializer_list<WGroupId> xs) {
    return std::unordered_set<WGroupId>(xs.begin(), xs.end());
}
static std::unordered_set<WGroupId> set_from(const std::vector<WGroupId>& xs) {
    return std::unordered_set<WGroupId>(xs.begin(), xs.end());
}

int main() {
    std::cout << "== DecayParent UNIT ==\n";

    ObsWilsonHelper(true);

    auto builder_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> builder = builder_spy;

    auto p_spy = std::make_shared<SpyObsParameterProxy>();
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_sm   = p_spy;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_flav = p_spy;

    auto qcd = std::make_shared<FakeObsQCDProxy>();
    auto use_marty_api = std::make_shared<FakeCoreBool>(false);

    auto freezer_spy = std::make_shared<SpyWilsonFreezer>();
    std::shared_ptr<IWilsonFreezer<WGroupId>> freezer = freezer_spy;

    ObservablePortsConfig ports(builder, p_sm, p_flav, qcd, use_marty_api, freezer);

    DummyDecay d(DecayMapper::to_id(Decays::B__l_l), 160.0, 4.8, QCDOrder::NONE, ports);

    {
        freezer_spy->clear();
        p_spy->calls_pid.clear();

        d.enable();

        assert(d.is_enabled());
        assert(builder_spy->build_calls == 1);
        assert(builder_spy->last_groups == set_from({GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime)}));

        assert(builder_spy->proxy->basis_called);
        assert(builder_spy->proxy->last_basis == WilsonBasis::B_STANDARD);

        assert(d.load_called);
        assert(!p_spy->calls_pid.empty());
    }

    {
        int prev_build = builder_spy->build_calls;
        d.load_called = false;
        builder_spy->proxy->basis_called = false;

        d.enable();

        assert(builder_spy->build_calls == prev_build);
        assert(!builder_spy->proxy->basis_called);
        assert(!d.load_called);
    }

    {
        d.disable();
        assert(!d.is_enabled());
    }

    {
        d.set_order(QCDOrder::NNLO);
        assert(d.current_order() == QCDOrder::NLO);

        d.set_order(QCDOrder::LO);
        assert(d.current_order() == QCDOrder::NLO);
    }

    {
        ObsWilsonHelper(true);
        use_marty_api->v = true;

        DummyDecay d2(DecayMapper::to_id(Decays::B__l_l), 160.0, 4.8, QCDOrder::NONE, ports);
        d2.set_order(QCDOrder::NLO);
        assert(d2.current_order() == QCDOrder::LO);
    }

    std::cout << "UNIT OK\n";
    return 0;
}
