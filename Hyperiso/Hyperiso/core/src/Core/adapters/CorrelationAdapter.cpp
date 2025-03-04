#include "CorrelationAdapter.h"

template CorrelationMatrixPair<ParamId> CorrelationAdapter::from_db_node<ParamId>(std::shared_ptr<Node>);
template CorrelationMatrixPair<Observables> CorrelationAdapter::from_db_node<Observables>(std::shared_ptr<Node>);

template<typename T>
CorrelationMatrixPair<T> CorrelationAdapter::from_db_node(std::shared_ptr<Node> root) {
    CorrelationMatrixPair<T> corr_matrices;
    for (auto &list_item : root->getGroup({"correlations"})) {
        emplace_correlation(corr_matrices, std::get<std::shared_ptr<Node>>(list_item.second));
    }
    return corr_matrices;
}

void CorrelationAdapter::emplace_correlation(CorrelationMatrixPair<ParamId> &corr_matrices, std::shared_ptr<Node> leaf) {
    auto parse_pid = [leaf] (size_t i) { 
        LhaID id {std::get<std::string>(leaf->get("id_" + std::to_string(i)))};
        auto block = std::get<std::string>(leaf->get("block_" + std::to_string(i)));
        return ParamId {ParamRouter::GetType(block, id), block, id};
    };

    if (!leaf->contains("block_1") || !leaf->contains("block_2") || !leaf->contains("id_1") 
        || !leaf->contains("id_2") || !leaf->contains("stat_correlation")) {
        LOG_ERROR("CorrelationAdapter", "Node doesn't have all necessary keys for parameter correlation.");
    }

    ParamId pid_1 = parse_pid(1);
    ParamId pid_2 = parse_pid(2);
    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = leaf->contains("syst_correlation") ? std::get<double>(leaf->get("syst_correlation")) : 0;

    corr_matrices.emplace(std::make_pair(pid_1, pid_2), stat_value, syst_value);
}

void CorrelationAdapter::emplace_correlation(CorrelationMatrixPair<Observables> &corr_matrices, std::shared_ptr<Node> leaf) {
    if (!leaf->contains("id_1") || !leaf->contains("id_2") || !leaf->contains("stat_correlation")) {
        LOG_ERROR("CorrelationAdapter", "Node doesn't have all necessary keys for observable correlation.");
    }

    Observables obs_1 = ObservableMapper::enum_elt(std::get<std::string>(leaf->get("id_1")));
    Observables obs_2 = ObservableMapper::enum_elt(std::get<std::string>(leaf->get("id_2")));
    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = leaf->contains("syst_correlation") ? std::get<double>(leaf->get("syst_correlation")) : 0;

    corr_matrices.emplace(std::make_pair(obs_1, obs_2), stat_value, syst_value);
}
