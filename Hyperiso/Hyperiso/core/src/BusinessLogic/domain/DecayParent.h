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

using std::chrono::high_resolution_clock;
using std::chrono::duration;

struct IFormFactorConfig {
    virtual ~IFormFactorConfig() = default;
};

class DecayParent {

protected:
    std::map<Observables, std::shared_ptr<OperatorNode>> roots;
    QCDOrder max_order = QCDOrder::LO; //DEFAULT AS LO, using default at least once, need to check (error if NONE)
    std::shared_ptr<ObsWilsonBuilder> w_builder;
    std::shared_ptr<ObsWilsonProxy> w_proxy;
    WilsonBuildConfig w_config {};
    bool enabled {false};

    QCDOrder check_max_order(QCDOrder order) const;

public:
    virtual ~DecayParent() = default;
    DecayParent(double matching_scale, double hadronic_scale, QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder);

    void enable();
    void disable();
    void set_order(QCDOrder new_order);

    scalar_t compute_observable(Observables obs);
    size_t get_n_evals(Observables obs);

    std::shared_ptr<OperatorNode> get_wilson_node(ScaleType scale=ScaleType::MATCHING, WilsonBasis basis=WilsonBasis::B_STANDARD);

    virtual void build_op_tree() = 0;

    template<typename EnumType>
    void set_config_flag(EnumType e) {
        throw std::logic_error("This flag is not supported for this decay");
    };
};

template <typename ConcreteDecay, typename... SupportedEnums>
class ConfigurableDecayParent : public DecayParent {
public:
    ConfigurableDecayParent(double matching_scale, double hadronic_scale,
                            QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder)
        : DecayParent(matching_scale, hadronic_scale, order, wilson_builder)
    {}

    template <typename EnumType>
    void set_config_flag(EnumType e) {
        if constexpr ((std::is_same_v<EnumType, SupportedEnums> || ...)) {
            static_cast<ConcreteDecay*>(this)->set_config_flag_c(e);
        } else {
            throw std::logic_error("Unsupported flag type for this decay");
        }
    }
};

#endif // __DECAYPARENT_H__