#include "DecayParent.h"
#include <iostream>
#include <algorithm>

QCDOrder DecayParent::check_max_order(QCDOrder order) const {
    if (order > max_order) {
        LOG_WARN("QCD order for decay cannot be higher than", OrderMapper::str(max_order));
        return max_order;
    }
    return order;
}

DecayParent::DecayParent(DecayId id, double matching_scale, double hadronic_scale, QCDOrder order, ObservablePortsConfig& ports) : 
    ports(ports), use_marty(ports.iobs_use_marty), p(ports.iobspp_sm), iobs_qcdp(ports.iobs_qcdp), iobs_wfreezer(ports.iobs_wfreezer), w_helper(ports.iobs_whelper) {
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
    if (this->enabled) {
        load_params();
        return;
    }
    if (!w_helper) {
        w_helper = std::make_shared<ObsWilsonHelper>();
    }
    w_helper->build(this->w_config, this->w_builder, iobs_wfreezer);
    this->w_proxy = this->w_builder->get_proxy();
    this->w_proxy->set_basis(WilsonBasis::B_STANDARD);
    load_params();
    this->enabled = true;
}

void DecayParent::disable() {
    this->enabled = false;
}

void DecayParent::set_order(QCDOrder new_order) {
    if (use_marty->get() && new_order > QCDOrder::LO) {
        LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
        new_order = QCDOrder::LO;
    }

    QCDOrder effective_order = check_max_order(new_order);
    if (this->w_config.order == effective_order) {
        return;
    }

    if (this->enabled) {
        LOG_WARN(
            "Changing QCD order for enabled decay",
            DecayMapper::str(this->id),
            "from",
            OrderMapper::str(this->w_config.order),
            "to",
            OrderMapper::str(effective_order),
            ": cached parameters are invalidated."
        );
        disable();
    }

    this->w_config.order = effective_order;
}

void DecayParent::set_n_threads(size_t n_threads) {
    LOG_WARN(
        "Thread configuration is not implemented for decay",
        DecayMapper::str(this->id),
        ". Requested",
        n_threads,
        "threads."
    );
}

size_t DecayParent::get_n_threads() const {
    return 1;
}

bool DecayParent::supports_thread_config() const {
    return false;
}

void DecayParent::set_bins(std::vector<std::pair<double, double>> new_bins) {
    for (auto& [qmin, qmax] : new_bins)
        if (qmin > qmax)
            LOG_WARN("Found inverted bin [", qmin, ",", qmax, "].");

    this->bins = new_bins;
}

void DecayParent::add_bin(std::pair<double, double> new_bin) {
    if (new_bin.first > new_bin.second)
        LOG_WARN("Found inverted bin [", new_bin.first, ",", new_bin.second, "].");

    if (!this->bins.has_value()) {
        this->bins = std::vector<std::pair<double, double>>{};
    }

    auto& current_bins = this->bins.value();
    if (std::find(current_bins.begin(), current_bins.end(), new_bin) == current_bins.end()) {
        current_bins.emplace_back(new_bin);
    }
}

std::optional<std::vector<std::pair<double, double>>> DecayParent::get_bins() {
    return this->bins;
}
