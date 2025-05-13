#include "ParamBlockLoader.h"

void ParamBlockLoader::load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) {
    LOG_INFO("Loading parameter blocks from", src_file.string());
    auto np = NodeProviderFactory::createNodeProvider(src_file);
    std::cout << "herhe" << std::endl;
    auto src = np->provide_db_as_node();
    // src->printJSON();
    std::cout << "ahah" << std::endl;
    for (auto &bk : src->get_keys()) {
        auto block = std::make_shared<Block>();
        block->blockname = bk;
        for (auto &vk : src->getGroup({bk})) {
            auto node = std::get<std::shared_ptr<Node>>(vk.second);

            if (!node->contains("central_value")) {
                LOG_ERROR("ParamBlockLoader", "Node doesn't have all necessary keys.");
            }

            auto value = node->get("central_value");
            if (std::holds_alternative<BlockName>(value)) {
                continue;
            } else {
                auto val = std::get<double>(value);
                auto stat = node->contains("stat_error") ? node->get("stat_error") : 0.;
                auto syst = node->contains("syst_error") ? node->get("syst_error") : 0.;
                block->store(LhaID(vk.first), std::make_shared<Parameter>(Parameter(ParamId(bk, std::string(vk.first)), val, std::get<double>(stat), std::get<double>(syst))));  //TODO : just be careful to string conversion, maybe better things to do.
            }
        }

        dest->emplace(bk, block);
    }

    LOG_INFO("Parameter blocks loaded");
}