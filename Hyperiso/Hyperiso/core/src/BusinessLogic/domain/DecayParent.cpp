#include "DecayParent.h"
#include <iostream>

QCDOrder DecayParent::check_max_order(QCDOrder order) const {
    if (order > max_order) {
        LOG_WARN("QCD order for decay cannot be higher than", OrderMapper::str(max_order));
        return max_order;
    }
    return order;
}

DecayParent::DecayParent(double matching_scale, double hadronic_scale, QCDOrder order, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) {
    this->w_config.matching_scale = matching_scale;
    this->w_config.hadronic_scale = hadronic_scale;
    this->w_config.order = check_max_order(order);
    this->w_builder = wilson_builder;
}

void DecayParent::enable() {
    // TODO : Manage enabling properly (don't call it each time eval() is called)
    // if (this->enabled) {
    //     return;
    // }

    ObsWilsonHelper::build(this->w_config, this->w_builder);
    this->w_proxy = this->w_builder->get_proxy();
    this->w_proxy->set_basis(WilsonBasis::B_STANDARD);
    build_op_tree();
    // this->enabled = true;
}

void DecayParent::disable() {
    this->enabled = false;
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
    auto t1 = high_resolution_clock::now();
    auto _ = truc->calculate();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms = t2 - t1;
    LOG_DEBUG("Evaluation of observable", ObservableMapper::str(obs), "took", ms.count(), "ms");
    return _;
}

size_t DecayParent::get_n_evals(Observables obs) {
    return roots.at(obs)->get_n_evals();
}

std::shared_ptr<OperatorNode> DecayParent::get_wilson_node(ScaleType scale, WilsonBasis basis) {
    auto wilson_node = std::make_shared<OperatorNode>("wilson", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 0; });

    for (WGroup group: this->w_config.groups) {
        for (WCoef c: WCoefMapper::get_group(group)) {
            std::string storage_block = GroupMapper::str(group, scale, basis);
            for (size_t order=1; order <= (size_t)this->w_config.order; order++) {
                LhaID c_id = WCoefMapper::flha_full(c, (QCDOrder)order, ContributionType::TOTAL);
                auto c_node = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, storage_block, c_id));
                wilson_node->addChild(c_node);
            }
        }
    }

    return wilson_node;
}
