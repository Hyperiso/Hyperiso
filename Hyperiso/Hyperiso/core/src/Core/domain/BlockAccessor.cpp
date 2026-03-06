#include "BlockAccessor.h"

scalar_t BlockAccessor::getValue(const BlockName& blockName, LhaID id) const {
    return this->at(blockName)->retrieve(id)->get_val();
}

std::shared_ptr<Parameter> BlockAccessor::getParameter(const BlockName& blockName, LhaID id) const {
    return this->at(blockName)->retrieve(id);
}

bool BlockAccessor::has_param(const BlockName& blockName, LhaID id) const {
    if (!this->contains(blockName)) return false;
    return this->at(blockName)->contains(id);
}

void BlockAccessor::setValue(const BlockName& blockName, LhaID id, scalar_t value) {
    if (!this->contains(blockName)) {
        throw std::invalid_argument("Block not found " + blockName.to_string());
    }
    auto& blk = this->at(blockName);
    if (blk->contains(id)) {
        blk->assign(id, value);
    } else {
        blk->store(id, std::make_shared<Parameter>(ParamId(blockName, id), value, 0., 0.));
    }
}

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

std::unordered_set<BlockName> BlockAccessor::get_block_names() const {
    std::unordered_set<BlockName> out;
    out.reserve(key_to_name_.size());
    for (auto& [k, bn] : key_to_name_) out.insert(bn);
    return out;
}


void BlockAccessor::remove_item(const BlockName& block_name, LhaID id) {
    if (this->contains(block_name)) {
        this->at(block_name)->remove(id);
    } else {
        LOG_WARN("Cannot remove item from non-existing block", block_name);
    }
}


std::shared_ptr<BlockAccessor> BlockAccessor::operator[](std::unordered_set<BlockName> block_names) {
    auto sub_block_accessor = std::make_shared<BlockAccessor>();

    for (const auto& block_name : block_names) {
        if (!this->contains(block_name)) {
            for (auto elem : *this) {
                std::cout << elem.first << std::endl;
            }
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
    const auto& blk = this->at(block_name);
    if (!blk->has_scale()) {
        LOG_ERROR("Block", block_name, "has no scale");
    }
    return blk->get_scale();
}

bool BlockAccessor::has_scale(const BlockName& block_name) const {
    return this->at(block_name)->has_scale();
}

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

            auto [statL, systL] = pLhs->get_std();
            bool lhsHasNoUncert = (statL == 0 && systL == 0);

            if (!lhsHasNoUncert) {
                res->setParameter(b, id, pLhs);
                continue;
            }

            auto [statR, systR] = pRhs->get_std();

            auto merged = std::make_shared<Parameter>(*pLhs);
            merged->set_std(statR, systR); 

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

std::unordered_map<std::string, std::shared_ptr<Block>>
BlockAccessor::get_block_sources(const BlockName& block_name) const {
    return this->at(block_name)->get_source_blocks();
}

std::unordered_map<ParamId, std::shared_ptr<Parameter>>
BlockAccessor::get_parameter_sources(const BlockName& block_name, LhaID id) const {
    return this->at(block_name)->retrieve(id)->get_source_parameters();
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



void BlockAccessor::emplace(const BlockName& name, std::shared_ptr<Block> block) {
    std::string existing_key = key_for(name);

    std::string key;
    if (!existing_key.empty()) {
        key = existing_key;
    } else {
        std::string best;
        bool first = true;
        for (const auto& a : name.get_alias()) {
            auto an = normalize(a);
            if (first || an < best) { best = an; first = false; }
        }
        key = best;
    }

    base_t::operator[](key) = std::move(block);

    auto& stored = base_t::operator[](key);
    if (stored) stored->bind_self(stored);
    
    merge_name_into_key(key, name);
}

void BlockAccessor::detach_block(const BlockName& block_name) {
    auto blk = this->at(block_name);

    if (auto dep = std::dynamic_pointer_cast<DependentBlock>(blk)) {
        dep->detach();
    }
}

void BlockAccessor::detach_parameter(const BlockName& block_name, LhaID id) {
    auto p = this->at(block_name)->retrieve(id);

    if (auto dep = std::dynamic_pointer_cast<DependentParameter>(p)) {
        dep->detach();
    }
}

std::string BlockAccessor::normalize(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : std::string(s)) out.push_back(char(std::toupper(c)));
    return out;
}

std::string BlockAccessor::key_for(std::string_view alias) const {
    const std::string a = normalize(alias);

    if (auto it = alias_to_key_.find(a); it != alias_to_key_.end())
        return it->second;

    if (base_t::find(a) != base_t::end())
        return a;

    const std::string raw{alias};
    if (base_t::find(raw) != base_t::end())
        return raw;

    return "";
}

std::string BlockAccessor::key_for(const std::string& alias) const {
    auto a = normalize(alias);
    auto it = alias_to_key_.find(a);
    if (it == alias_to_key_.end()) return "";
    return it->second;
}

std::string BlockAccessor::key_for(const BlockName& name) const {
    for (const auto& a : name.get_alias()) {
        auto k = key_for(a);
        if (!k.empty()) return k;
    }
    auto one = name.to_string();
    if (!one.empty()) {
        auto k = key_for(one);
        if (!k.empty()) return k;
    }
    return "";
}

std::string BlockAccessor::resolve_key(std::string_view alias) const {
    return key_for(alias);
}

std::string BlockAccessor::resolve_key(const BlockName& name) const {
    return key_for(name);
}

std::string BlockAccessor::choose_key(const std::unordered_set<std::string>& aliases_norm) {
    if (aliases_norm.empty()) return "";
    return *std::min_element(aliases_norm.begin(), aliases_norm.end());
}



void BlockAccessor::merge_name_into_key(const std::string& key, const BlockName& name) {
    auto& full = key_to_name_[key]; // crée si absent (vide)
    for (const auto& a : name.get_alias()) full.addAlias(a);

    for (const auto& a : full.get_alias()) {
        alias_to_key_[ normalize(a) ] = key;
    }
}

bool BlockAccessor::contains(std::string_view name) const {
    return !resolve_key(name).empty();
}

bool BlockAccessor::contains(const BlockName& name) const {
    return !resolve_key(name).empty();
}

std::shared_ptr<Block>& BlockAccessor::at(const BlockName& name) {
    auto key = resolve_key(name);
    if (key.empty())
        throw std::invalid_argument("Block not found " + name.to_string());
    return std::unordered_map<std::string, std::shared_ptr<Block>>::at(key);
}

const std::shared_ptr<Block>& BlockAccessor::at(const BlockName& name) const {
    auto key = resolve_key(name);
    if (key.empty())
        throw std::invalid_argument("Block not found " + name.to_string());
    return std::unordered_map<std::string, std::shared_ptr<Block>>::at(key);
}

std::shared_ptr<Block>& BlockAccessor::at(std::string_view name) {
    auto key = resolve_key(name);
    if (key.empty())
        throw std::invalid_argument("Block not found " + std::string(name));
    return std::unordered_map<std::string, std::shared_ptr<Block>>::at(key);
}

const std::shared_ptr<Block>& BlockAccessor::at(std::string_view name) const {
    auto key = resolve_key(name);
    if (key.empty())
        throw std::invalid_argument("Block not found " + std::string(name));
    return std::unordered_map<std::string, std::shared_ptr<Block>>::at(key);
}

bool BlockAccessor::has_param(std::string_view block, const LhaID& id) const {
    if (!contains(block)) return false;
    return at(block)->contains(id);
}

scalar_t BlockAccessor::getValue(std::string_view block, const LhaID& id) const {
    return at(block)->retrieve(id)->get_val();
}

void BlockAccessor::setValue(std::string_view block, const LhaID& id, scalar_t value) {
    auto& blk = at(block);
    if (blk->contains(id)) {
        blk->assign(id, value);
    } else {
        blk->store(id, std::make_shared<Parameter>(ParamId(std::string(block), id), value, 0, 0));
    }
}

void BlockAccessor::remove_item(std::string_view block, const LhaID& id) {
    at(block)->remove(id);
}

bool BlockAccessor::has_scale(std::string_view block) const {
    if (!contains(block)) return false;
    return at(block)->has_scale();
}

double BlockAccessor::get_scale(std::string_view block) const {
    return at(block)->get_scale();
}

void BlockAccessor::erase_block(std::string_view alias_or_key) {
    const std::string key = key_for(alias_or_key);
    if (key.empty()) return;  

    std::unordered_map<std::string, std::shared_ptr<Block>>::erase(key);

    auto itn = key_to_name_.find(key);
    if (itn != key_to_name_.end()) {
        for (const auto& a : itn->second.get_alias()) {
            auto ita = alias_to_key_.find(normalize(a));
            if (ita != alias_to_key_.end() && ita->second == key) {
                alias_to_key_.erase(ita);
            }
        }
        key_to_name_.erase(itn);
    }
}

void BlockAccessor::erase_block(const BlockName& name) {
    for (const auto& a : name.get_alias()) {
        const std::string key = key_for(a);
        if (!key.empty()) {
            erase_block(std::string_view{key});
            return;
        }
    }
}