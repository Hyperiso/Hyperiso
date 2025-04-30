#include "DecayParent.h"
#include <iostream>

QCDOrder DecayParent::check_max_order(QCDOrder order) const {
    if (order > max_order) {
        LOG_WARN("QCD order for decay cannot be higher than", OrderMapper::str(max_order));
        return max_order;
    }
    return order;
}

DecayParent::DecayParent(double matching_scale, double hadronic_scale, QCDOrder order, std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> wilson_builder) {
    this->w_config.matching_scale = matching_scale;
    this->w_config.hadronic_scale = hadronic_scale;
    this->w_config.order = check_max_order(order);
    this->w_builder = wilson_builder;
}

void DecayParent::enable() {
    ObsWilsonHelper::build(this->w_config, this->w_builder);
    this->w_proxy = this->w_builder->get_proxy();
    build_op_tree();
}

//TODO : everything here, just for make it works
void DecayParent::set_order(QCDOrder new_order) {
    if (ObsUseMarty().get() && new_order > QCDOrder::LO) {
        LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
        new_order = QCDOrder::LO;
    }
    
    if (this->w_config.order == new_order) {
        return;
    }
    
    if (this->w_config.order == QCDOrder::NONE) {
        this->w_config.order = check_max_order(new_order);
    } else {
        LOG_WARN("QCD order for this decay has already been set.");
    }
}

scalar_t DecayParent::compute_observable(Observables obs) {
    // TODO
    auto truc = roots.at(obs);
    auto _ = truc->calculate();
    return _;
}

size_t DecayParent::get_n_evals(Observables obs) {
    return roots.at(obs)->get_n_evals();
}
