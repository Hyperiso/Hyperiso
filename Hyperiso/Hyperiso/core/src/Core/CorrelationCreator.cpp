#include "CorrelationCreator.h"

template<typename T>
CorrelationMatrixPair<T> CorrelationCreator::from_db_node(std::shared_ptr<Node> root) {
    CorrelationMatrixPair corr_matrices;
    for (auto &list_item : root->getGroup({"correlations"})) {
        emplace_correlation(corr_matrices, std::get<std::shared_ptr<Node>>(list_item.second));
    }
    return corr_matrices;
}

void CorrelationCreator::emplace_correlation(CorrelationMatrixPair<ParamId> &corr_matrices, std::shared_ptr<Node> leaf) {
    auto parse_pid = [leaf] (size_t i) { 
        LhaID id {std::get<std::string>(leaf->get("id_" + std::to_string(i)))};
        auto block = std::get<std::string>(leaf->get("block_" + std::to_string(i)));
        return ParamId {ParamRouter::GetType(block, id), block, id};
    };

    if (!leaf->contains("block_1") || !leaf->contains("block_2") || !leaf->contains("id_1") 
        || !leaf->contains("id_2") || !leaf->contains("stat_correlation") || !leaf->contains("syst_correlation")) {
        LOG_ERROR("CorrelationCreator", "Node doesn't have all necessary keys for parameter correlation.");
    }

    ParamId pid_1 = parse_pid(1);
    ParamId pid_2 = parse_pid(2);
    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = std::get<double>(leaf->get("syst_correlation"));

    corr_matrices.stat.emplace(std::make_pair(pid_1, pid_2), stat_value);
    corr_matrices.stat.emplace(std::make_pair(pid_2, pid_1), stat_value);
    corr_matrices.syst.emplace(std::make_pair(pid_1, pid_2), syst_value);
    corr_matrices.syst.emplace(std::make_pair(pid_2, pid_1), syst_value);
}

void CorrelationCreator::emplace_correlation(CorrelationMatrixPair<Observables> &corr_matrices, std::shared_ptr<Node> leaf) {
    if (!leaf->contains("id_1") || !leaf->contains("id_2") || !leaf->contains("stat_correlation") || !leaf->contains("syst_correlation")) {
        LOG_ERROR("CorrelationCreator", "Node doesn't have all necessary keys for observable correlation.");
    }

    Observables obs_1 = ObservableMapper::enum_elt(std::get<std::string>(leaf->get("id_1")));
    Observables obs_2 = ObservableMapper::enum_elt(std::get<std::string>(leaf->get("id_2")));
    auto value = std::get<double>(leaf->get("value"));

    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = std::get<double>(leaf->get("syst_correlation"));

    corr_matrices.stat.emplace(std::make_pair(obs_1, obs_2), stat_value);
    corr_matrices.stat.emplace(std::make_pair(obs_2, obs_1), stat_value);
    corr_matrices.syst.emplace(std::make_pair(obs_1, obs_2), syst_value);
    corr_matrices.syst.emplace(std::make_pair(obs_2, obs_1), syst_value);
}
