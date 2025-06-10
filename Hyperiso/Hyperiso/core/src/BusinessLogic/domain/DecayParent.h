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
    DecayParent(double matching_scale, double hadronic_scale, QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder);

    void enable();
    void disable();
    void set_order(QCDOrder new_order);

    scalar_t compute_observable(Observables obs);
    size_t get_n_evals(Observables obs);

    std::shared_ptr<OperatorNode> get_wilson_node();

    virtual void build_op_tree() = 0;
};

#endif // __DECAYPARENT_H__