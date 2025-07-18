#include "ParamBlockLoader.h"

void ParamBlockLoader::load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) {
    LOG_INFO("Loading parameter blocks from", src_file.string());
    auto np = NodeProviderFactory::createNodeProvider(src_file);
    auto src = np->provide_db_as_node();
    
    for (auto &bk : src->get_keys()) {
        auto block = std::make_shared<Block>();
        block->blockname = bk;
        LOG_INFO("Loading block", bk);
        for (auto &vk : src->getGroup({bk})) {
            auto node = std::get<std::shared_ptr<Node>>(vk.second);

            if (!node->contains("central_value")) {
                LOG_ERROR("ParamBlockLoader", "Node doesn't have all necessary keys.");
            }

            if (node->contains("scale")) {
                auto scale = node->get("scale");
                if (!block->has_scale()) {
                    block->set_scale(std::get<double>(scale));
                }
            }

            auto value = node->get("central_value");
            if (std::holds_alternative<BlockName>(value)) {
                continue;
            } else {
                auto val = std::get<double>(value);
                auto stat = node->contains("stat_error") ? node->get("stat_error") : 0.;
                auto syst = node->contains("syst_error") ? node->get("syst_error") : 0.;
                

                block->store(LhaID(vk.first), std::make_shared<Parameter>(Parameter(ParamId(bk, std::string(vk.first)), val, std::get<double>(stat), std::get<double>(syst)))); 
            }
        }

        dest->emplace(bk, block);
    }

    LOG_INFO("Parameter blocks loaded");
}

void ParamBlockLoader::save(fs::path dest_file, std::shared_ptr<BlockAccessor> src){
    LOG_INFO("Saving blocks to ", dest_file.string());

    auto db_manager = std::make_shared<DBManager>();
    auto nd = std::make_shared<Node>();
    auto block_names = src->get_block_names();

    for (const auto &block_name : block_names) {
        auto block = src->at(block_name);
        auto block_ids = getAllIDs();

        for (const auto &id : block_ids) {

        }

        db_manager->write_to_file(dest_file, nd);
    }
}