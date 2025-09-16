#include "CorrelationAdapter.h"

template void CorrelationLoader<ParamId>::load(std::shared_ptr<CorrelationMatrixPair<ParamId>>, fs::path);
template void CorrelationLoader<ObservableId>::load(std::shared_ptr<CorrelationMatrixPair<ObservableId>>, fs::path);

template <typename T>
void CorrelationLoader<T>::load(std::shared_ptr<CorrelationMatrixPair<T>> dest, fs::path src_file) { 
    auto np = NodeProviderFactory::createNodeProvider(src_file);
    auto src = np->provide_db_as_node();

    // src->printJSON();

    for (auto &list_item : src->getGroup({"correlations"})) {
        emplace_correlation(dest, std::get<std::shared_ptr<Node>>(list_item.second));
    }
}

template<>
void CorrelationLoader<ParamId>::emplace_correlation(std::shared_ptr<CorrelationMatrixPair<ParamId>> corr_matrices, std::shared_ptr<Node> leaf) {
    auto parse_pid = [leaf] (size_t i) { 
        Node::Value id_value = leaf->get("id_" + std::to_string(i));
        LhaID id;
        if (std::holds_alternative<int>(id_value)) {
            id = LhaID{std::get<int>(id_value)};
        } else {
            id = LhaID{std::get<BlockName>(id_value)};
        }
        auto block = std::get<BlockName>(leaf->get("block_" + std::to_string(i)));
        return ParamId {ParamRouter::GetType(block, id), block, id};
    };

    if (!leaf->contains("block_1") || !leaf->contains("block_2") || !leaf->contains("id_1") 
        || !leaf->contains("id_2") || !leaf->contains("stat_correlation")) {
        LOG_ERROR("CorrelationLoader", "Node doesn't have all necessary keys for parameter correlation.");
    }

    ParamId pid_1 = parse_pid(1);
    ParamId pid_2 = parse_pid(2);
    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = leaf->contains("syst_correlation") ? std::get<double>(leaf->get("syst_correlation")) : 0;

    corr_matrices->emplace(std::make_pair(pid_1, pid_2), stat_value, syst_value);
}

template<>
void CorrelationLoader<ObservableId>::emplace_correlation(std::shared_ptr<CorrelationMatrixPair<ObservableId>> corr_matrices, std::shared_ptr<Node> leaf) {
    if (!leaf->contains("id_1") || !leaf->contains("id_2") || !leaf->contains("stat_correlation")) {
        LOG_ERROR("CorrelationLoader", "Node doesn't have all necessary keys for observable correlation.");
    }
    //TODO : do some checking with value()
    ObservableId obs_1 = ObservableMapper::from_flha(LhaID(std::get<BlockName>(leaf->get("id_1")))).value();
    ObservableId obs_2 = ObservableMapper::from_flha(LhaID(std::get<BlockName>(leaf->get("id_2")))).value();
    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = leaf->contains("syst_correlation") ? std::get<double>(leaf->get("syst_correlation")) : 0;

    corr_matrices->emplace(std::make_pair(obs_1, obs_2), stat_value, syst_value);
}
