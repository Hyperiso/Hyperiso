#include "DecayParent.h"
#include <iostream>

QCDOrder DecayParent::check_max_order(QCDOrder order) const {
    if (order > max_order) {
        LOG_WARN("QCD order for decay cannot be higher than", OrderMapper::str(max_order));
        return max_order;
    }
    return order;
}

DecayParent::DecayParent(double matching_scale, double hadronic_scale, QCDOrder order) {
    this->w_config.matching_scale = matching_scale;
    this->w_config.hadronic_scale = hadronic_scale;
    this->w_config.order = check_max_order(order);
}

void DecayParent::enable() {
    ObsWilsonHelper::build(this->w_config);
    this->w_proxy = std::make_shared<ObsWilsonProxy>();
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
    std::cout << ObservableMapper::str(obs) << std::endl;

    std::cout << "----------------------------" << std::endl;
    for (auto& elem : roots) {
        std::cout << ObservableMapper::str(elem.first) << std::endl;
    }
    auto truc = roots.at(obs);
    std::cout << "its was fine" << std::endl;
    auto _ = truc->calculate();
    std::cout << "its fine" << std::endl;
    return _;
}

size_t DecayParent::get_n_evals(Observables obs) {
    return roots.at(obs)->get_n_evals();
}
