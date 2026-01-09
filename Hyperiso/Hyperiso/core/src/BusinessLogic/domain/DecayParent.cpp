#include "DecayParent.h"
#include <iostream>

QCDOrder DecayParent::check_max_order(QCDOrder order) const {
    if (order > max_order) {
        LOG_WARN("QCD order for decay cannot be higher than", OrderMapper::str(max_order));
        return max_order;
    }
    return order;
}

DecayParent::DecayParent(DecayId id, double matching_scale, double hadronic_scale, QCDOrder order, ObservablePortsConfig& ports) : 
    ports(ports), use_marty(ports.iobs_use_marty), p(ports.iobspp_sm), iobs_qcdp(ports.iobs_qcdp), iobs_wfreezer(ports.iobs_wfreezer) {
    bind_wilson_builder(ports.iobswb);
    this->id = id;
    this->w_config.matching_scale = matching_scale;
    this->w_config.hadronic_scale = hadronic_scale;
    this->w_config.order = check_max_order(order);
}

void DecayParent::bind_wilson_builder(std::shared_ptr<IObsWilsonBuilder> &wilson_builder) {
    this->w_builder = wilson_builder;
}

void DecayParent::enable() {
    if (this->enabled)
        return;
    ObsWilsonHelper::build(this->w_config, this->w_builder, iobs_wfreezer);
    this->w_proxy = this->w_builder->get_proxy();
    this->w_proxy->set_basis(WilsonBasis::B_STANDARD);
    load_params();
    this->enabled = true;
}

void DecayParent::disable() {
    this->enabled = false;
}

//TODO : everything here, just for make it works
void DecayParent::set_order(QCDOrder new_order) {
    if (use_marty->get() && new_order > QCDOrder::LO) {
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
