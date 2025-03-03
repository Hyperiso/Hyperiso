#include "BlockAccessor.h"

void BlockAccessor::addBlock(const std::string &name,
                             std::shared_ptr<Block> block) {
    blocks[name] = block;
}

double BlockAccessor::getValue(const std::string& blockName, LhaID id) const {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getValue(id);
    }
    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(id));
}

Parameter BlockAccessor::getParameter(const std::string &blockName, LhaID id) const {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getParameter(id);
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

void BlockAccessor::setParameter(const std::string &blockName, LhaID id, const Parameter &source) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        it->second->setParameter(id, source);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::map<LhaID, double> BlockAccessor::getAllValues(std::string blockName)
{
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getAllValues();
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::unordered_set<std::string> BlockAccessor::get_block_names() {
    return get_keys(blocks);
}

std::shared_ptr<Block> BlockAccessor::get_block(const std::string &block_name) {
    if (this->blocks.contains(block_name)) {
        return this->blocks[block_name];
    }
    LOG_ERROR("BlockAccessor", "Cannot retrieve non-existing block", block_name);
}

void BlockAccessor::remove_block(const std::string &block_name) {
    if (this->blocks.contains(block_name)) {
        this->blocks.erase(block_name);
    } else {
        LOG_WARN("Cannot remove non-existing block", block_name);
    }
}

void BlockAccessor::remove_item(const std::string &block_name, LhaID id) {
    if (this->blocks.contains(block_name)) {
        this->blocks[block_name]->remove_parameter(id);
    } else {
        LOG_WARN("Cannot remove item from non-existing block", block_name);
    }
}

std::shared_ptr<BlockAccessor> BlockAccessor::operator[](std::unordered_set<std::string> block_names) {
    auto sub_block_accessor = std::make_shared<BlockAccessor>();

    for (const auto& block_name : block_names) {
        if (!this->has_block(block_name)) {
            LOG_ERROR("BlockAccessor", "Block", block_name, "doesn't exist. Cannot extract.");
        }
        sub_block_accessor->addBlock(block_name, this->get_block(block_name));
    }

    return sub_block_accessor;
}

std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::make_shared<BlockAccessor>();
    for (const auto &b : lhs->get_block_names()) {
        res->addBlock(b, std::make_shared<MapBlock>(lhs->get_block(b)));
    }

    for (const auto &b : rhs->get_block_names()) {
        if (res->has_block(b)) {
            LOG_ERROR("BlockAccessor", "Cannot merge blocks with common blocks using no-priority operator +. Use >> for priority-merge.");
        }
        res->addBlock(b, std::make_shared<MapBlock>(rhs->get_block(b)));
    }

    return res;
}

std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::make_shared<BlockAccessor>();
    for (const auto &b : rhs->get_block_names()) {
        res->addBlock(b, std::make_shared<MapBlock>(rhs->get_block(b)));
    }

    for (const auto &b : lhs->get_block_names()) {
        if (res->has_block(b)) {
            for (const auto& id : lhs->get_block(b)->getAllIDs()) {
                res->setValue(b, id, lhs->getValue(b, id), true);
            }
        } else {
            res->addBlock(b, std::make_shared<MapBlock>(lhs->get_block(b)));
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

void BlockAccessor::addDependentBlock(const std::string& name, std::shared_ptr<DependentBlock>& dependant_block, const std::string& sourceName, std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateFunc) {
    auto sourceBlock = get_block(sourceName);
    if (!sourceBlock) {
        throw std::invalid_argument("Source block not found");
    }

    dependant_block = std::make_shared<DependentBlock>(sourceBlock, recalculateFunc);
    dependant_block->blockname = name;
    dependant_block->init();
    blocks[name] = dependant_block;
}