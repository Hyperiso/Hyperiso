#include "ParamBlockLoader.h"

// void ParamBlockLoader::load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) {
//     LOG_INFO("Loading parameter blocks from", src_file.string());
//     auto np = NodeProviderFactory::createNodeProvider(src_file);
//     auto src = np->provide_db_as_node();
    
//     for (auto &bk : src->get_keys()) {
//         auto block = std::make_shared<Block>();
//         block->blockname = bk;
//         LOG_INFO("Loading block", bk);
//         for (auto &vk : src->getGroup({bk})) {
//             auto node = std::get<std::shared_ptr<Node>>(vk.second);

//             if (!node->contains("central_value")) {
//                 LOG_ERROR("ParamBlockLoader", "Node doesn't have all necessary keys.");
//             }

//             if (node->contains("scale")) {
//                 auto scale = node->get("scale");
//                 if (!block->has_scale()) {
//                     block->set_scale(std::get<double>(scale));
//                 }
//             }

//             auto value = node->get("central_value");
//             if (std::holds_alternative<BlockName>(value)) {
//                 continue;
//             } else {
//                 auto val = std::get<double>(value);
//                 auto stat = node->contains("stat_error") ? node->get("stat_error") : 0.;
//                 auto syst = node->contains("syst_error") ? node->get("syst_error") : 0.;
                

//                 block->store(LhaID(vk.first), std::make_shared<Parameter>(Parameter(ParamId(bk, std::string(vk.first)), val, std::get<double>(stat), std::get<double>(syst)))); 
//             }
//         }

//         dest->emplace(bk, block);
//     }

//     LOG_INFO("Parameter blocks loaded");
// }

void ParamBlockLoader::load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) {
    LOG_INFO("Loading parameter blocks from", src_file.string());
    auto np  = NodeProviderFactory::createNodeProvider(src_file);
    auto src = np->provide_db_as_node();

    for (auto &bk : src->get_keys()) {
        auto block = std::make_shared<Block>();
        block->blockname = bk;
        LOG_INFO("Loading block", bk);

        auto group = src->getGroup({bk});  // map<BlockName, Node::Value>

        for (auto &vk : group) {
            const auto& key = vk.first;
            const auto& val = vk.second;

            // 1) gérer la scale du bloc si présente ici
            if (key == "scale") {
                if (!block->has_scale()) {
                    if (std::holds_alternative<double>(val)) {
                        block->set_scale(std::get<double>(val));
                    } else if (std::holds_alternative<int>(val)) {
                        block->set_scale(static_cast<double>(std::get<int>(val)));
                    } else {
                        LOG_WARN("ParamBlockLoader", "Non-numeric block scale under ", bk);
                    }
                }
                continue; // très important: ne pas traiter "scale" comme une entrée
            }

            // 2) les vraies entrées doivent être des noeuds
            if (!std::holds_alternative<std::shared_ptr<Node>>(val)) {
                LOG_WARN("ParamBlockLoader", "Unexpected non-node entry under ", bk, " key ", key, " — skipping");
                continue;
            }

            auto node = std::get<std::shared_ptr<Node>>(val);

            if (!node->contains("central_value")) {
                LOG_ERROR("ParamBlockLoader", "Node doesn't have all necessary keys (central_value).");
                continue;
            }

            // 3) si une scale par-entrée existe, on peut l'ignorer ou la traiter selon ton design
            if (node->contains("scale") && !block->has_scale()) {
                auto scale = node->get("scale");
                if (std::holds_alternative<double>(scale))
                    block->set_scale(std::get<double>(scale));
                else if (std::holds_alternative<int>(scale))
                    block->set_scale(static_cast<double>(std::get<int>(scale)));
            }

            // 4) extraire la valeur centrale et erreurs
            auto value = node->get("central_value");
            if (std::holds_alternative<BlockName>(value)) {
                continue; // ignore les labels éventuels
            }

            double val_central = std::get<double>(value);

            auto stat = node->contains("stat_error") ? node->get("stat_error") : Node::Value{0.0};
            auto syst = node->contains("syst_error") ? node->get("syst_error") : Node::Value{0.0};

            double stat_d = std::holds_alternative<double>(stat) ? std::get<double>(stat)
                             : std::holds_alternative<int>(stat) ? static_cast<double>(std::get<int>(stat))
                                                                 : 0.0;
            double syst_d = std::holds_alternative<double>(syst) ? std::get<double>(syst)
                             : std::holds_alternative<int>(syst) ? static_cast<double>(std::get<int>(syst))
                                                                 : 0.0;

            block->store(
                LhaID(vk.first),
                std::make_shared<Parameter>(Parameter(ParamId(bk, std::string(vk.first)),
                                                      val_central, stat_d, syst_d))
            );
        }

        dest->emplace(bk, block);
    }

    LOG_INFO("Parameter blocks loaded");
}
