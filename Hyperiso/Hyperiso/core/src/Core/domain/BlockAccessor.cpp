#include "BlockAccessor.h"

double BlockAccessor::getValue(const std::string& blockName, LhaID id) const {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return it->second->getValue(id);
    }
    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(id));
}

Parameter BlockAccessor::getParameter(const std::string &blockName, LhaID id) const {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return it->second->getParameter(id);
    }
    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(id));
}

bool BlockAccessor::exist(const std::string blockName, LhaID id) const {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return true;
    }
    throw false;
}

void BlockAccessor::setValue(const std::string& blockName, LhaID id, double value, bool force) {
    auto it = this->find(blockName);
    if (it != this->end()) {
        it->second->setValue(id, value, force);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

void BlockAccessor::setMode(const std::string& blockName, LhaID id, ParameterMode mode) {
    auto it = this->find(blockName);
    if (it != this->end()) {
        it->second->setMode(id, mode);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

void BlockAccessor::setParameter(const std::string &blockName, LhaID id, const Parameter &source) {
    auto it = this->find(blockName);
    if (it != this->end()) {
        it->second->setParameter(id, source);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::map<LhaID, double> BlockAccessor::getAllValues(std::string blockName) {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return it->second->getAllValues();
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::unordered_set<std::string> BlockAccessor::get_block_names() {
    return get_keys(*this);
}

void BlockAccessor::remove_item(const std::string &block_name, LhaID id) {
    if (this->contains(block_name)) {
        this->at(block_name)->remove_parameter(id);
    } else {
        LOG_WARN("Cannot remove item from non-existing block", block_name);
    }
}

std::shared_ptr<BlockAccessor> BlockAccessor::operator[](std::unordered_set<std::string> block_names) {
    auto sub_block_accessor = std::make_shared<BlockAccessor>();

    for (const auto& block_name : block_names) {
        if (!this->contains(block_name)) {
            LOG_ERROR("BlockAccessor", "Block", block_name, "doesn't exist. Cannot extract.");
        }
        sub_block_accessor->emplace(block_name, this->at(block_name));
    }

    return sub_block_accessor;
}

std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::make_shared<BlockAccessor>();
    for (const auto &b : lhs->get_block_names()) {
        res->emplace(b, std::make_shared<MapBlock>(lhs->at(b)));
    }

    for (const auto &b : rhs->get_block_names()) {
        if (res->contains(b)) {
            LOG_ERROR("BlockAccessor", "Cannot merge blocks with common blocks using no-priority operator +. Use >> for priority-merge.");
        }
        res->emplace(b, std::make_shared<MapBlock>(rhs->at(b)));
    }

    return res;
}

std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::make_shared<BlockAccessor>();
    for (const auto &b : rhs->get_block_names()) {
        res->emplace(b, std::make_shared<MapBlock>(rhs->at(b)));
    }

    for (const auto &b : lhs->get_block_names()) {
        if (res->contains(b)) {
            for (const auto& id : lhs->at(b)->getAllIDs()) {
                res->setValue(b, id, lhs->getValue(b, id), true);
            }
        } else {
            res->emplace(b, std::make_shared<MapBlock>(lhs->at(b)));
        }
    }

    return res;
}

std::ostream &operator<<(std::ostream &os, std::shared_ptr<BlockAccessor> ba) {
    for (auto& k : ba->get_block_names()) {
        os << "Block " << k << ":\n";
        for (auto &[id, val] : ba->at(k)->getAllValues()) {
            os << '\t' << id << ": " << val << '\n';
        }
        os << '\n';
    }
    return os;
}

void BlockAccessor::addDependentBlock(const std::string& name, std::shared_ptr<DependentBlock>& dependant_block, const std::string& sourceName, std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateFunc) {
    auto sourceBlock = at(sourceName);
    if (!sourceBlock) {
        throw std::invalid_argument("Source block not found");
    }

    dependant_block = std::make_shared<DependentBlock>(sourceBlock, recalculateFunc);
    dependant_block->blockname = name;
    dependant_block->init();
    this->emplace(name, dependant_block);
}