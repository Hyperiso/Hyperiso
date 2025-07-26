#include "ParamBlockLoader.h"
#include "IDataLoader.h"
template class IDataLoader<BlockAccessor>;

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
    auto node = std::make_shared<Node>();

    for (const auto &block_name : src->get_block_names()) {
        const auto &block = src->at(block_name);
        const auto &items = block->getItems();

        std::map<BlockName, Node::Value> block_data;
        for (const auto &[id, param] : items) {
            auto node_param = std::make_shared<Node>();
                node_param->set(param->get_val(), "central_value");
            // if (param->has_stat_error()) {
            //     node_param->set(param->get_stat_error(), "stat_error");
            // }
            // if (param->has_syst_error()) {
            //     node_param->set(param->get_syst_error(), "syst_error");
            // }
            if (block->has_scale()) {
                node_param->set(block->get_scale(), "scale");
            }

            block_data[id.to_string()] = node_param;
        }
        node->setGroup({block_name}, block_data);
        LOG_INFO("Saving block_name", block_name);
        LOG_INFO("Block size", block_data.size());
        //node->printJSON(2);
    }
    //LOG_INFO("Saving node", node->get_keys());

    db_manager->write_to_file(dest_file, node);
}