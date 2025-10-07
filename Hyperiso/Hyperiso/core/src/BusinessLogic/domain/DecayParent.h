#ifndef DECAYPARENT_H
#define DECAYPARENT_H

#include <map>
#include <string>
#include "General.h"
#include "WilsonInterface.h"
#include "Node.h"
#include "ObsUseMarty.h"
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

using std::chrono::high_resolution_clock;
using std::chrono::duration;

class DecayParent  {

protected:
    QCDOrder max_order = QCDOrder::LO;
    std::shared_ptr<ObsWilsonBuilder> w_builder;
    std::shared_ptr<ObsWilsonProxy> w_proxy;
    WilsonBuildConfig w_config {};
    bool enabled {false};
    DecayId id;

    QCDOrder check_max_order(QCDOrder order) const;

public:
    virtual ~DecayParent() = default;
    DecayParent(DecayId custom_id, double matching_scale, double hadronic_scale, QCDOrder order);
    DecayParent(DecayId custom_id, double matching_scale, double hadronic_scale, QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder);
    
    void bind_wilson_builder(std::shared_ptr<ObsWilsonBuilder>& wilson_builder);
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
                            QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder)
        : DecayParent(id, matching_scale, hadronic_scale, order, wilson_builder)
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
                            QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder)
        : DecayParent(id, matching_scale, hadronic_scale, order, wilson_builder)
    {}

    DecayParentConfigurable(DecayId id, double matching_scale, double hadronic_scale,
                            QCDOrder order)
        : DecayParent(id, matching_scale, hadronic_scale, order)
    {}

    void set_config(std::any) final override {

    }
};

#endif // __DECAYPARENT_H__