#include "CorrelationAdapter.h"
#include "ExperimentObs.h"

template void CorrelationLoader<ParamId>::load(
    std::shared_ptr<CorrelationMatrixPair<ParamId>>, fs::path, bool);

// template void CorrelationLoader<ExperimentObs>::load(
//     std::shared_ptr<CorrelationMatrixPair<ExperimentObs>>, fs::path, bool);

template <typename T>
void CorrelationLoader<T>::load(std::shared_ptr<CorrelationMatrixPair<T>> dest,
                                fs::path src_file,
                                bool block_in_blocks) {
    auto np = DBNodeProviderFactory::createDBNodeProvider(src_file);
    auto src = np->provide_db_as_node();

    for (auto &list_item : src->getGroup({"correlations"})) {
        emplace_correlation(dest, std::get<std::shared_ptr<DBNode>>(list_item.second));
    }
}

template<>
void CorrelationLoader<ParamId>::emplace_correlation(
    std::shared_ptr<CorrelationMatrixPair<ParamId>> corr_matrices,
    std::shared_ptr<DBNode> leaf) {

    auto parse_pid = [leaf](size_t i) {
        DBNode::Value id_value = leaf->get("id_" + std::to_string(i));
        LhaID id;

        if (std::holds_alternative<int>(id_value)) {
            id = LhaID{std::get<int>(id_value)};
        } else {
            id = LhaID{std::get<BlockName>(id_value)};
        }

        auto block = std::get<BlockName>(leaf->get("block_" + std::to_string(i)));
        return ParamId{ParamRouter::GetType(block, id), block, id};
    };

    if (!leaf->contains("block_1") || !leaf->contains("block_2")
        || !leaf->contains("id_1") || !leaf->contains("id_2")
        || !leaf->contains("stat_correlation")) {
        LOG_ERROR("CorrelationLoader",
                  "DBNode doesn't have all necessary keys for parameter correlation.");
        return;
    }

    ParamId pid_1 = parse_pid(1);
    ParamId pid_2 = parse_pid(2);

    auto stat_value = std::get<double>(leaf->get("stat_correlation"));
    auto syst_value = leaf->contains("syst_correlation")
                    ? std::get<double>(leaf->get("syst_correlation"))
                    : 0.0;

    corr_matrices->emplace(std::make_pair(pid_1, pid_2), stat_value, syst_value);
}

template<>
void CorrelationLoader<ExperimentObs>::load(
    std::shared_ptr<CorrelationMatrixPair<ExperimentObs>> dest,
    fs::path src_file,
    bool block_in_blocks) {

    auto np = DBNodeProviderFactory::createDBNodeProvider(src_file);
    auto src = np->provide_db_as_node();

    if (!src->contains("correlations")) {
        LOG_ERROR("CorrelationLoader",
                  "DBNode doesn't contain 'correlations' for observable correlations.");
        return;
    }

    // Nouveau format :
    // {
    //   "correlations": {
    //      "DEFAULT":  [ ... ],
    //      "DEFAULT2": [ ... ]
    //   }
    // }

    for (auto &exp_item : src->getGroup({"correlations"})) {
        const std::string experiment = exp_item.first;

        for (auto &list_item : src->getGroup({"correlations", experiment})) {
            auto leaf = std::get<std::shared_ptr<DBNode>>(list_item.second);

            if (!leaf->contains("id_1") || !leaf->contains("id_2")
                || !leaf->contains("stat_correlation")) {
                LOG_ERROR("CorrelationLoader",
                          "DBNode doesn't have all necessary keys for observable correlation.");
                continue;
            }

            BinnedObservableId obs_1 =
                BinnedObservableId::from_flha(LhaID(std::get<BlockName>(leaf->get("id_1"))));
            BinnedObservableId obs_2 =
                BinnedObservableId::from_flha(LhaID(std::get<BlockName>(leaf->get("id_2"))));

            auto stat_value = std::get<double>(leaf->get("stat_correlation"));
            auto syst_value = leaf->contains("syst_correlation")
                            ? std::get<double>(leaf->get("syst_correlation"))
                            : 0.0;

            dest->emplace(
                std::make_pair(
                    ExperimentObs{experiment, obs_1},
                    ExperimentObs{experiment, obs_2}
                ),
                stat_value,
                syst_value
            );
        }
    }
}