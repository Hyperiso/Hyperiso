#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <any>

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
#include "IObsParameterProxy.h"
#include "IObsWilsonBuilder.h"
#include "IObsWilsonProxy.h"
#include "IObsCoreAPI.h"
#include "IObsQCDProxy.h"
#include "IWilsonFreezer.h"
#include "DecayParent.h"
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

class DummyFreezer : public IWilsonFreezer<WGroupId> {
public:
    void freeze(WGroupId) override {}
    void unfreeze(WGroupId) override {}
};

class DummyObsWilsonProxy : public IObsWilsonProxy {
public:
    complex_t getM (WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getFM(WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getR (WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }
    complex_t getFR(WGroup, WCoef, QCDOrder, ContributionType) override { return {0,0}; }

    std::map<QCDOrder, complex_t> getSM(WGroup, WCoef, ContributionType) override { return {}; }
    std::map<QCDOrder, complex_t> getSR(WGroup, WCoef, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAM (WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAR (WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAFM(WGroup, QCDOrder, ContributionType) override { return {}; }
    std::map<WCoef, complex_t> getAFR(WGroup, QCDOrder, ContributionType) override { return {}; }

    std::shared_ptr<IObsWilsonBuilder> get_builder() override { return nullptr; }
    std::unordered_set<WilsonBasis> get_bases(WGroupId) override { return {WilsonBasis::B_STANDARD}; }
    void set_basis(WilsonBasis) override {}
};

class DummyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    std::shared_ptr<DummyObsWilsonProxy> proxy = std::make_shared<DummyObsWilsonProxy>();
    void build(std::shared_ptr<AbstractConfig>) override {}
    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return proxy; }
};

class SpyParamProxy : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    std::unordered_map<ParamId, scalar_t> values;

    mutable int calls = 0;
    mutable std::vector<ParamId> called_pids;

    scalar_t operator()(const ParamId& pid, DataType) override {
        calls++;
        called_pids.push_back(pid);
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

class SimpleDecay : public DecayParent {
public:
    SimpleDecay(DecayId did, ObservablePortsConfig& ports)
      : DecayParent(did, 0.0, 0.0, QCDOrder::NONE, ports) {}
    void load_params() override {}
    std::vector<ObservableValue> compute_observable(Observables) override { return {}; }
    std::vector<ObservableValue> compute_observable(ObservableId) override { return {}; }
    void set_config(std::any) override {}
};

class BindCheckDecay : public DecayParent {
public:
    explicit BindCheckDecay(DecayId did, ObservablePortsConfig& ports)
        : DecayParent(did, 0.0, 0.0, QCDOrder::LO, ports) {}

    bool is_bound_to(const std::shared_ptr<IObsWilsonBuilder>& b) const {
        return w_builder == b;
    }

    void load_params() override {}
    std::vector<ObservableValue> compute_observable(Observables) override { return {}; }
    std::vector<ObservableValue> compute_observable(ObservableId) override { return {}; }
    void set_config(std::any) override {}
};

static inline scalar_t R(double x){ return scalar_t{x,0.0}; }

static bool contains_pid(const std::vector<ParamId>& v, const ParamId& pid) {
    for (auto& x : v) if (x == pid) return true;
    return false;
}

int main() {
    std::cout << "== ObsManager UNIT ==\n";

    auto wb_impl = std::make_shared<DummyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> wb = wb_impl; 

    auto pp_spy = std::make_shared<SpyParamProxy>();
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_sm   = pp_spy;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_flav = pp_spy;

    auto qcd = std::make_shared<FakeObsQCDProxy>();
    auto use_marty = std::make_shared<FakeCoreBool>(false);
    auto freezer = std::make_shared<DummyFreezer>();

    ParamId pid_mW(ParameterType::SM, "MASS", 24);
    ParamId pid_mb(ParameterType::SM, "QCD", LhaID(5,2));
    pp_spy->values[pid_mW] = R(80.0);
    pp_spy->values[pid_mb] = R(4.2);

    ObservablePortsConfig ports(wb, p_sm, p_flav, qcd, use_marty, freezer);

    pp_spy->calls = 0;
    pp_spy->called_pids.clear();

    ObsManager mgr(ports, false);

    auto did = DecayMapper::to_id(Decays::B__l_l);
    mgr.add_custom_decay(did, std::make_shared<SimpleDecay>(did, ports));

    assert(pp_spy->calls >= 2);
    assert(contains_pid(pp_spy->called_pids, pid_mW));
    assert(contains_pid(pp_spy->called_pids, pid_mb));

    {
        const auto oid = BinnedObservableId(ObservableMapper::to_id(Observables::BR_BS_MUMU));

        mgr.add_obs(Observables::BR_BS_MUMU, QCDOrder::LO, /*add_deps*/false);

        auto cur = mgr.get_current_obss();
        assert(std::count(cur.begin(), cur.end(), oid) == 1);

        auto obs_ptr = mgr.get_obs(Observables::BR_BS_MUMU);
        assert(obs_ptr != nullptr);
        assert(obs_ptr->getId() == oid);
    }

    {
        const auto oid = BinnedObservableId(ObservableMapper::to_id(Observables::BR_BS_MUMU));
        mgr.remove_obs(Observables::BR_BS_MUMU);

        auto cur = mgr.get_current_obss();
        assert(std::count(cur.begin(), cur.end(), oid) == 1);
    }

    {
        auto did = DecayMapper::to_id(Decays::B__l_l);

        auto custom = std::make_shared<BindCheckDecay>(did, ports);

        mgr.add_custom_decay(did, custom);

        assert(custom->is_bound_to(wb));
    }

    std::cout << "UNIT OK\n";
    return 0;
}
