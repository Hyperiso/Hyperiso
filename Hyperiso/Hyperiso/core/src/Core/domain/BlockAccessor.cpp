#include "BlockAccessor.h"

scalar_t BlockAccessor::getValue(const BlockName& blockName, LhaID id) const {
    if (this->contains(blockName)) {
        return this->at(blockName)->retrieve(id)->get_val();
    }

    throw std::invalid_argument("Block " + blockName + " not found (LhaID: " + id.to_string() + ")");
}

std::shared_ptr<Parameter> BlockAccessor::getParameter(const BlockName& blockName, LhaID id) const {
    if (this->contains(blockName)) {
        return this->at(blockName)->retrieve(id);
    }

    throw std::invalid_argument("Block " + blockName + " not found (LhaID: " + id.to_string() + ").");
}

bool BlockAccessor::has_param(const BlockName& blockName, LhaID id) const {
    return this->contains(blockName) && this->at(blockName)->contains(id);
}

void BlockAccessor::setValue(const BlockName& blockName, LhaID id, scalar_t value) {
    if (this->contains(blockName)) {
        if (this->at(blockName)->contains(id)) {
            this->at(blockName)->assign(id, value);
        } else {
            this->at(blockName)->store(id, std::make_shared<Parameter>(ParamId(blockName, id), value, 0., 0.));
        }
    } else {
        throw std::invalid_argument("Block not found " + blockName);
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
    if (this->contains(blockName)) {
        this->at(blockName)->store_or_assign(id, source);
    } else {
        throw std::invalid_argument("Block not found : " + blockName);
    }
}

std::map<LhaID, scalar_t> BlockAccessor::getAllValues(BlockName blockName) {
    if (this->contains(blockName)) {
        std::map<LhaID, scalar_t> values;
        for(auto& [id, p] : this->at(blockName)->getItems()) {
            values.emplace(id, p->get_val());
        }
        return values;
    } else {
        throw std::invalid_argument("Block not found : " + blockName);
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

double BlockAccessor::get_scale(const BlockName& block_name) const {
    auto it = std::find_if(
        this->begin(),
        this->end(),
        [&](const auto& pair) { return pair.first == block_name; }
    );

    if (it == this->end()) {
        std::cout << std::make_shared<BlockAccessor>(*this) << std::endl;
        LOG_ERROR("Block", block_name, "not found in BlockAccessor");
    }

    if (it->second->has_scale()) {
        return it->second->get_scale();
    } else {
        LOG_ERROR("Block", block_name, "has no scale");
    }
    
}

bool BlockAccessor::has_scale(const BlockName& block_name) const {
    auto it = std::find_if(
        this->begin(),
        this->end(),
        [&](const auto& pair) { return pair.first == block_name; }
    );

    if (it == this->end()) {
        std::cout << std::make_shared<BlockAccessor>(*this) << std::endl;
        LOG_ERROR("Block", block_name, "not found in BlockAccessor");
    }

    return it->second->has_scale();
}


// std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs) {
//     auto res = std::make_shared<BlockAccessor>();
//     for (const auto &b : rhs->get_block_names()) {
//         res->emplace(b, std::make_shared<Block>(rhs->at(b)));
//     }

//     for (const auto &b : lhs->get_block_names()) {
//         if (res->contains(b)) {
//             for (const auto& id : lhs->at(b)->getAllIDs()) {
//                 res->setParameter(b, id, lhs->getParameter(b, id));
//             }
//         } else {
//             res->emplace(b, std::make_shared<Block>(lhs->at(b)));
//         }
//     }

//     return res;
// }

std::shared_ptr<BlockAccessor> operator>>(
    std::shared_ptr<BlockAccessor> lhs,
    std::shared_ptr<BlockAccessor> rhs
) {
    auto res = std::make_shared<BlockAccessor>();

    for (const auto& b : rhs->get_block_names()) {
        res->emplace(b, std::make_shared<Block>(rhs->at(b)));
    }

    for (const auto& b : lhs->get_block_names()) {
        auto lhsBlock = lhs->at(b);

        if (!res->contains(b)) {
            res->emplace(b, std::make_shared<Block>(*lhsBlock));
            continue;
        }

        auto rhsBlock = rhs->at(b);
        auto resBlock = res->at(b);

        const auto rhsIds = rhsBlock->getAllIDs();

        for (const auto& id : lhsBlock->getAllIDs()) {
            auto pLhs = lhs->getParameter(b, id);

            bool inRhs = std::find(rhsIds.begin(), rhsIds.end(), id) != rhsIds.end();

            if (!inRhs) {
                res->setParameter(b, id, pLhs);
                continue;
            }

            auto pRhs = rhs->getParameter(b, id);

            auto [statR, systR] = pRhs->get_std();
            bool rhsHasNoUncert = (statR == 0 && systR == 0);

            if (!rhsHasNoUncert) {
                continue;
            }

            auto [statL, systL] = pLhs->get_std();

            auto merged = std::make_shared<Parameter>(*pRhs);
            merged->set_std(statL, systL);

            res->setParameter(b, id, merged);
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

std::unordered_map<std::string, std::shared_ptr<Block>> BlockAccessor::get_block_sources(const BlockName& block_name) const {
    auto it = std::find_if(
        this->begin(),
        this->end(),
        [&](const auto& pair) { return pair.first == block_name; }
    );

    if (it == this->end()) {
        std::cout << std::make_shared<BlockAccessor>(*this) << std::endl;
        LOG_ERROR("Block", block_name, "not found in BlockAccessor");
    }
    return it->second->get_source_blocks();

}

std::unordered_map<ParamId, std::shared_ptr<Parameter>> BlockAccessor::get_parameter_sources(const BlockName& block_name, LhaID id) const {
    auto it = std::find_if(
        this->begin(),
        this->end(),
        [&](const auto& pair) { return pair.first == block_name; }
    );

    if (it == this->end()) {
        std::cout << std::make_shared<BlockAccessor>(*this) << std::endl;
        LOG_ERROR("Block", block_name, "not found in BlockAccessor");
    }
    return it->second->retrieve(id)->get_source_parameters();
}

std::unordered_set<ParamId>
BlockAccessor::get_all_source_parameters(const std::unordered_set<ParamId>& param_ids) const
{
    std::unordered_set<ParamId> result;
    std::unordered_set<ParamId> visited;

    std::function<void(const ParamId&)> dfs =
        [&](const ParamId& pid)
        {
            if (!visited.insert(pid).second) {
                return;
            }

            if (!this->has_param(pid.block, pid.code)) {
                result.insert(pid);
                return;
            }

            bool has_sources = false;

            auto param_sources = this->get_parameter_sources(pid.block, pid.code);
            if (!param_sources.empty()) {
                has_sources = true;
                for (const auto& [src_id, _] : param_sources) {
                    dfs(src_id);
                }
            }

            auto block_sources = this->get_block_sources(pid.block);
            if (!block_sources.empty()) {
                has_sources = true;

                for (const auto& [src_block_name, src_block] : block_sources) {
                    if (!src_block) continue;


                    const auto& items = src_block->getItems();
                    for (const auto& [lha_id, param_ptr] : items) {
                        if (!param_ptr) continue;
                        dfs(param_ptr->get_id());
                    }
                }
            }


            if (!has_sources) {
                result.insert(pid);
            }
        };

    for (const auto& pid : param_ids) {
        dfs(pid);
    }

    return result;
}