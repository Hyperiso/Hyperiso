#include "BlocksCreator.h"

std::shared_ptr<BlockAccessor> BlocksCreator::from_lha_reader(std::shared_ptr<LhaReader> reader) {
    auto ba = std::make_shared<BlockAccessor>();
    for (auto &bk : reader->getBlocksNames()) {
        auto block = std::make_shared<MapBlock>();
        block->blockname = bk;
        for (auto &vk : *reader->getBlock(bk)->getEntries()) {
            block->setValue(vk->getId(), reader->getValue<double>(bk, vk->getId()));
        }
    }

    return ba;
}

std::shared_ptr<BlockAccessor> BlocksCreator::from_db_node(std::shared_ptr<Node> root) {
    auto ba = std::make_shared<BlockAccessor>();
    for (auto &bk : root->get_keys()) {
        auto block = std::make_shared<MapBlock>();
        block->blockname = bk;
        for (auto &vk : root->getGroup({bk})) {
            auto value = std::get<std::shared_ptr<Node>>(vk.second)->get("central_value");
            block->setValue(std::stol(vk.first), std::get<double>(value));
            auto stat = std::get<std::shared_ptr<Node>>(vk.second)->get("stat_error");
            auto syst = std::get<std::shared_ptr<Node>>(vk.second)->get("syst_error");
            block->setDeviation(std::stol(vk.first), std::get<double>(stat), std::get<double>(syst));
        }
    }

    return ba;
}