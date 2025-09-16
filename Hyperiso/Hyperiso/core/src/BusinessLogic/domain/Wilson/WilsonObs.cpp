#include "WilsonObs.h"

complex_t WilsonDecay::get_matching_coef(WGroup group_id, WCoef coef_id) {
    return w_proxy->getFM(group_id, coef_id, w_config.order);
}

complex_t WilsonDecay::get_running_coef(WGroup group_id, WCoef coef_id) {
    return w_proxy->getFR(group_id, coef_id, w_config.order);
}

void WilsonDecay::build_op_tree() {
    auto wilson = this->get_wilson_node();

    auto C7 = std::make_shared<OperatorNode>("C7", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return get_matching_coef(WGroup::B, WCoef::C7); });
    C7->addChild(wilson);

    DecayMapper::register_custom("Wilson");
    ObservableMapper::register_custom("C7", {}, LhaID(42,42), "Wilson");
    ObservableId c7_id = ObservableMapper::id_of("C7");
    roots.emplace(std::pair{c7_id, C7});
}