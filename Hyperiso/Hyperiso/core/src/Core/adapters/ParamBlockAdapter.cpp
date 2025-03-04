#include "ParamBlockAdapter.h"

std::shared_ptr<BlockAccessor> ParamBlockAdapter::from_db_node(std::shared_ptr<Node> root) {
    auto ba = std::make_shared<BlockAccessor>();
    for (auto &bk : root->get_keys()) {
        auto block = std::make_shared<MapBlock>();
        block->blockname = bk;
        for (auto &vk : root->getGroup({bk})) {
            auto node = std::get<std::shared_ptr<Node>>(vk.second);
            if (!node->contains("central_value")) {
                LOG_ERROR("ParamBlockAdapter", "Node doesn't have all necessary keys.");
            }

            auto value = node->get("central_value");
            block->setValue(LhaID(vk.first), std::get<double>(value));

            auto stat = node->contains("stat_error") ? node->get("stat_error") : 0;
            auto syst = node->contains("syst_error") ? node->get("syst_error") : 0;
            block->setDeviation(LhaID(vk.first), std::get<double>(stat), std::get<double>(syst));
        }

        ba->addBlock(bk, block);
    }

    return ba;
}