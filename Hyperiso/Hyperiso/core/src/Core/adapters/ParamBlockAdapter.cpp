#include "ParamBlockAdapter.h"

void ParamBlockLoader::load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) {
    auto np = NodeProviderFactory::createNodeProvider(src_file);
    auto src = np->provide_db_as_node();
    
    for (auto &bk : src->get_keys()) {
        auto block = std::make_shared<MapBlock>();
        block->blockname = bk;
        for (auto &vk : src->getGroup({bk})) {
            auto node = std::get<std::shared_ptr<Node>>(vk.second);

            if (!node->contains("central_value")) {
                LOG_ERROR("ParamBlockLoader", "Node doesn't have all necessary keys.");
            }

            auto value = node->get("central_value");
            if (std::holds_alternative<std::string>(value)) {
                continue;
            } else {
                block->setValue(LhaID(vk.first), std::get<double>(value));
                auto stat = node->contains("stat_error") ? node->get("stat_error") : 0.;
                auto syst = node->contains("syst_error") ? node->get("syst_error") : 0.;
                block->setDeviation(LhaID(vk.first), std::get<double>(stat), std::get<double>(syst));  
            }
        }

        dest->addBlock(bk, block);
    }
}