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

#include "Observable.h"
#include "DecayParent.h"
#include "ObsWilsonHelper.h"
#include "ObsPortsConfig.h"
#include "IObsWilsonBuilder.h"
#include "IObsWilsonProxy.h"
#include "IObsParameterProxy.h"
#include "IObsCoreAPI.h"
#include "IObsQCDProxy.h"
#include "IWilsonFreezer.h"
#include "Configs.h"
#include "Parameter.h"
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

class SpyWilsonFreezer : public IWilsonFreezer<WGroupId> {
public:
    std::vector<WGroupId> froze, unfroze;
    void freeze(WGroupId g) override { froze.push_back(g); }
    void unfreeze(WGroupId g) override { unfroze.push_back(g); }
};

class SpyObsParamProxyObs : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    mutable int calls = 0;
    mutable ParamId last_pid{};
    mutable DataType last_dt = DataType::VALUE;

    std::map<std::pair<ParamId, DataType>, scalar_t> table;

    scalar_t operator()(const ParamId& pid, DataType dt) override {
        calls++;
        last_pid = pid;
        last_dt  = dt;
        auto it = table.find({pid, dt});
        if (it == table.end()) return scalar_t{0.0, 0.0};
        return it->second;
    }
    scalar_t operator()(const std::string&, const LhaID&, DataType) const override {
        return scalar_t{0.0, 0.0};
    }
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const override {
        return std::make_shared<Parameter>(pid, 0.0, 0.0, 0.0);
    }
};

class DummyObsParamProxy : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    scalar_t operator()(const ParamId&, DataType) override { return scalar_t{0.0, 0.0}; }
    scalar_t operator()(const std::string&, const LhaID&, DataType) const override { return scalar_t{0.0, 0.0}; }
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const override {
        return std::make_shared<Parameter>(pid, 0.0, 0.0, 0.0);
    }
};

class SpyObsWilsonProxy : public IObsWilsonProxy {
public:
    bool basis_called = false;
    WilsonBasis last_basis = WilsonBasis::B_STANDARD;

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

    void set_basis(WilsonBasis b) override {
        basis_called = true;
        last_basis = b;
    }
};

class SpyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    int build_calls = 0;
    std::shared_ptr<SpyObsWilsonProxy> proxy = std::make_shared<SpyObsWilsonProxy>();

    void build(std::shared_ptr<AbstractConfig>) override { build_calls++; }
    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return proxy; }
};

class TestDecay : public DecayParent {
public:
    mutable int load_calls = 0;
    mutable int compute_calls = 0;
    mutable ObservableId last_obs{};

    std::vector<ObservableValue> ret;

    TestDecay(DecayId did, ObservablePortsConfig& ports)
        : DecayParent(did, 160.0, 4.8, QCDOrder::LO, ports)
    {
        w_config.groups = { GroupMapper::to_id(WGroup::B) };
    }

    void load_params() override { load_calls++; }

    std::vector<ObservableValue> compute_observable(Observables) override { return {}; }

    std::vector<ObservableValue> compute_observable(ObservableId obs) override {
        compute_calls++;
        last_obs = obs;
        return ret;
    }

    void set_config(std::any) override {}
};

int main() {
    std::cout << "== Observable UNIT ==\n";

    ObsWilsonHelper(true);

    auto wb_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> wb = wb_spy;

    auto pp_dummy = std::make_shared<DummyObsParamProxy>();
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_sm = pp_dummy;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p_flav = pp_dummy;

    auto qcd = std::make_shared<FakeObsQCDProxy>();
    auto use_marty = std::make_shared<FakeCoreBool>(false);

    auto freezer_spy = std::make_shared<SpyWilsonFreezer>();
    std::shared_ptr<IWilsonFreezer<WGroupId>> freezer = freezer_spy;

    ObservablePortsConfig ports(wb, p_sm, p_flav, qcd, use_marty, freezer);

    auto decay = std::make_shared<TestDecay>(DecayMapper::to_id(Decays::B__l_l), ports);

    auto obspp_spy = std::make_shared<SpyObsParamProxyObs>();

    const ObservableId oid = ObservableMapper::to_id(Observables::BR_BS_MUMU);
    Observable obs(oid, decay, obspp_spy);

    {
        obspp_spy->calls = 0;

        const auto flha = ObservableMapper::binned_flha_of({oid, {0.,0.}}).value();
        ParamId pid(ParameterType::OBSERVABLE, "FOBS_DEFAULT", flha);

        obspp_spy->table[{pid, DataType::VALUE}] = scalar_t{1.23, 0.0};

        auto v = obs.get_exp_val();
        assert(std::abs(v.real() - 1.23) < 1e-12);
        assert(std::abs(v.imag()) < 1e-12);

        assert(obspp_spy->calls == 1);
        assert(obspp_spy->last_pid == pid);
        assert(obspp_spy->last_dt == DataType::VALUE);
    }

    {
        obspp_spy->calls = 0;

        const auto flha = ObservableMapper::binned_flha_of({oid, {0.,0.}}).value();
        ParamId pid(ParameterType::OBSERVABLE, "FOBS_DEFAULT", flha);
        const auto dt = UncertaintyTypeMapper::d_type(UncertaintyType::COMBINED);

        obspp_spy->table[{pid, dt}] = scalar_t{0.11, 0.0};

        auto u = obs.get_exp_uncertainty({0.0, 0.}, "DEFAULT", UncertaintyType::COMBINED);
        assert(std::abs(u.real() - 0.11) < 1e-12);
        assert(std::abs(u.imag()) < 1e-12);

        assert(obspp_spy->calls == 1);
        assert(obspp_spy->last_pid == pid);
        assert(obspp_spy->last_dt == dt);
    }

    {
        wb_spy->build_calls = 0;
        decay->load_calls = 0;
        decay->compute_calls = 0;
        wb_spy->proxy->basis_called = false;

        decay->ret = { ObservableValue(oid, 9.99) };

        auto out = obs.compute();

        // Current Observable semantics:
        // compute() is only a thin forwarding wrapper. It does not enable the decay.
        assert(wb_spy->build_calls == 0);
        assert(!wb_spy->proxy->basis_called);
        assert(decay->load_calls == 0);

        assert(decay->compute_calls == 1);
        assert(decay->last_obs == oid);

        assert(out.size() == 1);
        assert(out[0].id == oid);
        assert(std::abs(out[0].value - 9.99) < 1e-12);
    }

    {
        ParamId p1{ParameterType::SM, "MASS", 13};
        ParamId p2{ParameterType::SM, "SMINPUTS", 2};
        ParamId p3{ParameterType::FLAVOR, "FMASS", 531};

        obs.add_dependence(p1);
        obs.add_dependence(p2);
        obs.add_dependence(p3);

        assert(obs.get_dependences().count(p1) == 1);
        assert(obs.get_dependences().count(p2) == 1);
        assert(obs.get_dependences().count(p3) == 1);

        obs.add_dependence(p2);
        assert(obs.get_dependences().count(p2) == 1);
    }

    std::cout << "UNIT OK\n";
    return 0;
}
