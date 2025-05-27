#include "BlockAccessor.h"

scalar_t BlockAccessor::getValue(const BlockName& blockName, LhaID id) const {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return it->second->retrieve(id)->get_val();
    }

    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + id.to_string());
}

std::shared_ptr<Parameter> BlockAccessor::getParameter(const BlockName& blockName, LhaID id) const {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return it->second->retrieve(id);
    }

    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + id.to_string());
}

bool BlockAccessor::has_param(const BlockName& blockName, LhaID id) const {
    auto it = this->find(blockName);
    if (it != this->end()) {
        return it->second->contains(id);
    }
    return false;
}

void BlockAccessor::setValue(const BlockName& blockName, LhaID id, scalar_t value) {
    LOG_INFO("BlockAccessor::setBlockValue");
    auto it = this->find(blockName);
    if (it != this->end()) {
        if (it->second->contains(id)) {
            it->second->assign(id, value);
        } else {
            it->second->store(id, std::make_shared<Parameter>(ParamId(blockName, id), value, 0., 0.));
        }
    } else {
        throw std::invalid_argument("Block not found");
    }
}

// void BlockAccessor::setMode(const std::string& blockName, LhaID id, ParameterMode mode) {
//     auto it = this->find(blockName);
//     if (it != this->end()) {
//         it->second->setMode(id, mode);
//     } else {
//         throw std::invalid_argument("Block not found");
//     }
// }

void BlockAccessor::setParameter(const BlockName &blockName, LhaID id, std::shared_ptr<Parameter> source) {
    auto it = this->find(blockName);
    if (it != this->end()) {
        it->second->store(id, source);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::map<LhaID, double> BlockAccessor::getAllValues(BlockName blockName) {
    auto it = this->find(blockName);
    if (it != this->end()) {
        std::map<LhaID, double> values;
        for(auto& [id, p] : it->second->getItems()) {
            values.emplace(id, p->get_val());
        }
        return values;
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::unordered_set<BlockName> BlockAccessor::get_block_names() {
    return get_keys(*this);
}

void BlockAccessor::remove_item(const BlockName& block_name, LhaID id) {
    if (this->contains(block_name)) {
        this->at(block_name)->remove(id);
    } else {
        LOG_WARN("Cannot remove item from non-existing block", block_name);
    }
}

bool BlockAccessor::contains(const BlockName& block_name) const {
    auto keys = get_keys(*this);
    return std::any_of(
        keys.begin(),
        keys.end(),
        [&](const auto& bn) { return bn == block_name; }
    );
}

std::shared_ptr<Block> &BlockAccessor::at(const BlockName &block_name) {
    auto it = std::find_if(
        this->begin(),
        this->end(),
        [&](const auto& pair) { return pair.first == block_name; }
    );

    if (it == this->end()) {
        std::cout << std::make_shared<BlockAccessor>(*this) << std::endl;
        LOG_ERROR("Block", block_name, "not found in BlockAccessor");
    }

    return it->second;
}

const std::shared_ptr<Block> &BlockAccessor::at(const BlockName &block_name) const {
    auto it = std::find_if(
        this->begin(),
        this->end(),
        [&](const auto& pair) { return pair.first == block_name; }
    );

    if (it == this->end()) {
        std::cout << std::make_shared<BlockAccessor>(*this) << std::endl;
        LOG_ERROR("Block", block_name, "not found in BlockAccessor");
    }

    return it->second;
}

std::shared_ptr<BlockAccessor> BlockAccessor::operator[](std::unordered_set<BlockName> block_names) {
    auto sub_block_accessor = std::make_shared<BlockAccessor>();

    for (const auto& block_name : block_names) {
        if (!this->contains(block_name)) {
            LOG_ERROR("BlockAccessor", "Block", block_name, "doesn't exist. Cannot extract.");
        }
        auto block_ptr = this->at(block_name);
        sub_block_accessor->emplace(block_ptr->get_name(), block_ptr);
    }

    return sub_block_accessor;
}

std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::make_shared<BlockAccessor>();
    for (const auto &b : lhs->get_block_names()) {
        res->emplace(b, std::make_shared<Block>(lhs->at(b)));
    }

    for (const auto &b : rhs->get_block_names()) {
        if (res->contains(b)) {
            LOG_ERROR("BlockAccessor", "Cannot merge blocks with common blocks using no-priority operator +. Use >> for priority-merge.");
        }
        res->emplace(b, std::make_shared<Block>(rhs->at(b)));
    }

    return res;
}

std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
    auto res = std::make_shared<BlockAccessor>();
    for (const auto &b : rhs->get_block_names()) {
        res->emplace(b, std::make_shared<Block>(rhs->at(b)));
    }

    for (const auto &b : lhs->get_block_names()) {
        if (res->contains(b)) {
            for (const auto& id : lhs->at(b)->getAllIDs()) {
                res->setValue(b, id, lhs->getValue(b, id));
            }
        } else {
            res->emplace(b, std::make_shared<Block>(lhs->at(b)));
        }
    }

    return res;
}

std::ostream &operator<<(std::ostream &os, std::shared_ptr<BlockAccessor> ba) {
    for (auto& block_name : ba->get_block_names()) {
        os << "Block " << block_name << ":\n";
        for (auto &[id, val] : ba->getAllValues(block_name)) {
            os << '\t' << id << ": " << val << '\n';
        }
        os << '\n';
    }
    return os;
}