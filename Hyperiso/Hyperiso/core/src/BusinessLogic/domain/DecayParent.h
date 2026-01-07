#ifndef DECAYPARENT_H
#define DECAYPARENT_H

#include <map>
#include <string>
#include "Include.h"
#include "WilsonInterface.h"
#include "ObsWilsonBuilder.h"
#include "ObsWilsonProxy.h"
#include "ObsWilsonHelper.h"
#include "Math.h"
#include "Configs.h"
#include <chrono>
#include <type_traits>
#include <any>
#include "DefaultConfig.h"
#include "ObservableValue.h"
#include "ObsParameterProxy.h"
#include "ObsPortsConfig.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration;

class DecayParent  {

protected:
    QCDOrder max_order = QCDOrder::LO;
    ObservablePortsConfig& ports;
    std::shared_ptr<IObsWilsonBuilder> w_builder;
    std::shared_ptr<IObsWilsonProxy> w_proxy;
    std::shared_ptr<IObsCoreAPI<bool>> use_marty;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p;
    std::shared_ptr<IObsQCDProxy> iobs_qcdp;
    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer;
    WilsonBuildConfig w_config {};
    bool enabled {false};
    DecayId id;

    QCDOrder check_max_order(QCDOrder order) const;

public:
    virtual ~DecayParent() = default;
    // DecayParent(DecayId custom_id, double matching_scale, double hadronic_scale, QCDOrder order);
    DecayParent(DecayId custom_id, double matching_scale, double hadronic_scale, QCDOrder order, ObservablePortsConfig& ports);
    
    void bind_wilson_builder(std::shared_ptr<IObsWilsonBuilder>& wilson_builder);
    void enable();
    void disable();
    void set_order(QCDOrder new_order);

    virtual void load_params() = 0; 
    virtual std::vector<ObservableValue> compute_observable(Observables obs) = 0;
    virtual std::vector<ObservableValue> compute_observable(ObservableId obs) = 0;

    virtual void set_config(std::any cfg) = 0;

    DecayId get_id() { return id; };
};

template<typename T>
class DecayParentConfigurable : public DecayParent {
public:
    DecayParentConfigurable(DecayId id, double matching_scale, double hadronic_scale,
                            QCDOrder order, ObservablePortsConfig& ports)
        : DecayParent(id, matching_scale, hadronic_scale, order, ports)
    {}

     void set_config(std::any cfg) final override {
        if (auto p = std::any_cast<T>(&cfg)) { set_config_spe(*p); return; }
        if (auto p = std::any_cast<std::add_const_t<T>>(&cfg)) { set_config_spe(*p); return; }
        throw std::bad_any_cast{};
    }
    virtual void set_config_spe(T config) = 0;
};

template<>
struct DecayParentConfigurable<DecayConfig> : DecayParent {
public:
    DecayParentConfigurable(DecayId id, double matching_scale, double hadronic_scale,
                            QCDOrder order, ObservablePortsConfig& ports)
        : DecayParent(id, matching_scale, hadronic_scale, order, ports)
    {}

    // DecayParentConfigurable(DecayId id, double matching_scale, double hadronic_scale,
    //                         QCDOrder order)
    //     : DecayParent(id, matching_scale, hadronic_scale, order)
    // {}

    void set_config(std::any) final override {

    }
};

#endif // DECAYPARENT_H