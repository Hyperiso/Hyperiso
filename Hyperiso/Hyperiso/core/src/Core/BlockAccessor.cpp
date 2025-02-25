#include "BlockAccessor.h"

void BlockAccessor::addBlock(const std::string &name,
                             std::shared_ptr<Block> block) {
    blocks[name] = std::move(block);
}

double BlockAccessor::getValue(const std::string& blockName, LhaID id) const {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getValue(id);
    }
    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(id));
}

bool BlockAccessor::exist(const std::string blockName, LhaID id) const {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return true;
    }
    throw false;
}

bool BlockAccessor::has_block(const std::string blockName) const {
    return this->blocks.contains(blockName);
}

void BlockAccessor::setValue(const std::string& blockName, LhaID id, double value, bool force) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        it->second->setValue(id, value, force);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

void BlockAccessor::setMode(const std::string& blockName, LhaID id, ParameterMode mode) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        it->second->setMode(id, mode);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::map<LhaID, double> BlockAccessor::getAllValues(std::string blockName) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getAllValues();
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::vector<std::string> BlockAccessor::get_block_names() {
    std::vector<std::string> keys;
    for (auto key : blocks){
        keys.push_back(key.first);
    }
    return keys;
}

std::shared_ptr<Block> BlockAccessor::get_block(const std::string &block_name) {
    if (this->blocks.contains(block_name)) {
        return this->blocks[block_name];
    }
    LOG_ERROR("BlockAccessor", "Cannot retrieve non-existing block", block_name);
}

std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::shared_ptr<BlockAccessor>();
    for (const auto &b : lhs->get_block_names()) {
        res->addBlock(b, lhs->get_block(b));
    }

    for (const auto &b : rhs->get_block_names()) {
        if (res->has_block(b)) {
            LOG_ERROR("BlockAccessor", "Cannot merge blocks with common blocks using no-priority operator +. Use >> for priority-merge.");
        }
        res->addBlock(b, rhs->get_block(b));
    }

    return res;
}

std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::shared_ptr<BlockAccessor>();
    for (const auto &b : rhs->get_block_names()) {
        res->addBlock(b, rhs->get_block(b));
    }

    for (const auto &b : lhs->get_block_names()) {
        if (res->has_block(b)) {
            for (const auto& id : lhs->get_block(b)->getAllIDs()) {
                res->setValue(b, id, lhs->getValue(b, id), true);
            }
        } else {
            res->addBlock(b, lhs->get_block(b));
        }
    }

    return res;
}

std::ostream &operator<<(std::ostream &os, std::shared_ptr<BlockAccessor> ba) {
    for (auto& k : ba->get_block_names()) {
        os << "Block " << k << ":\n";
        for (auto &[id, val] : ba->get_block(k)->getAllValues()) {
            os << '\t' << id << ": " << val << '\n';
        }
        os << '\n';
    }
    return os;
}
