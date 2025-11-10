#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <memory>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"

struct BAKeyHash {
    std::size_t operator()(const std::pair<std::string,std::string>& k) const noexcept {
        std::size_t h1 = std::hash<std::string>()(k.first);
        std::size_t h2 = std::hash<std::string>()(k.second);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1<<6) + (h1>>2));
    }
};

class ParamOptimizer {
public:
    explicit ParamOptimizer(std::shared_ptr<BlockAccessor> scope)
        : scopes_{std::move(scope)} {}
    explicit ParamOptimizer(std::vector<std::shared_ptr<BlockAccessor>> scopes)
        : scopes_{std::move(scopes)} {}

    void set_value(const BlockName& block, const LhaID& id, scalar_t v) {
        ops_.push_back(OpSetValue{block, id, v});
    }
    void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) {
        ops_.push_back(OpSetParam{block, id, std::move(p)});
    }
    void remove(const BlockName& block, const LhaID& id) {
        ops_.push_back(OpRemove{block, id});
    }

    void commit(bool coalesce = true) {
        if (ops_.empty()) return;

        freeze_all_();

        const auto plan = coalesce ? coalesce_ops_() : ops_;

        // initial existence before modif (by block/id)
        std::unordered_map<std::pair<std::string,std::string>, bool, BAKeyHash> existed;
        existed.reserve(plan.size());
        for (const auto& v : plan) {
            std::visit([&](auto&& op){
                auto blk = find_block_(op.block);
                existed[{op.block, op.id.to_string()}] = blk->contains(op.id);
            }, v);
        }

        std::unordered_set<std::shared_ptr<Block>> blocks_needing_notify;

        for (const auto& v : plan) {
            std::visit([&](auto&& op) {
                using T = std::decay_t<decltype(op)>;
                auto blk = find_block_(op.block);
                const bool had = existed.at({op.block, op.id.to_string()});

                if constexpr (std::is_same_v<T, OpSetValue>) {
                    if (had) {
                        blk->assign(op.id, op.value);          // notify but frozen
                    } else {
                        blk->store(op.id, std::make_shared<Parameter>(ParamId(blk->get_name(), op.id), op.value, 0., 0.));
                        blocks_needing_notify.insert(blk);      // store does not notify, notify block at the end
                    }
                } else if constexpr (std::is_same_v<T, OpSetParam>) {
                    if (had) {
                        blk->assign(op.id, op.param);           // notify but frozen
                    } else {
                        blk->store(op.id, op.param);
                        blocks_needing_notify.insert(blk);
                    }
                } else if constexpr (std::is_same_v<T, OpRemove>) {
                    if (!had) {
                        //remove of key which wasn't here initialy, ignore
                        return;
                    }
                    // real remove, destroy obs, no notif planned then
                    blk->remove(op.id);
                }
            }, v);
        }

        //1 notif per block
        for (const auto& blk : blocks_needing_notify) {
            if (blk) blk->notifyObservers();
        }

        unfreeze_all_();
        ops_.clear();
    }

    void clear() { ops_.clear(); }

private:
    struct OpSetValue { BlockName block; LhaID id; scalar_t value; };
    struct OpSetParam { BlockName block; LhaID id; std::shared_ptr<Parameter> param; };
    struct OpRemove   { BlockName block; LhaID id; };
    using Op = std::variant<OpSetValue, OpSetParam, OpRemove>;

    std::vector<Op> coalesce_ops_() const {
        std::unordered_map<std::pair<std::string,std::string>, std::size_t, BAKeyHash> last;
        last.reserve(ops_.size());
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            const auto& v = ops_[i];
            std::visit([&](auto&& op){
                last[{op.block, op.id.to_string()}] = i; // last winning
            }, v);
        }
        std::vector<std::size_t> idx; idx.reserve(last.size());
        for (auto& kv : last) idx.push_back(kv.second);
        std::sort(idx.begin(), idx.end());
        std::vector<Op> res; res.reserve(idx.size());
        for (auto i : idx) res.push_back(ops_[i]);
        return res;
    }

    std::shared_ptr<Block> find_block_(const BlockName& name) const {
        std::shared_ptr<Block> found;
        bool ambiguous = false;
        for (const auto& ba : scopes_) {
            if (!ba) continue;
            if (ba->contains(name)) {
                auto b = ba->at(name);
                if (!found) found = b;
                else if (found.get() != b.get()) ambiguous = true;
            }
        }
        if (!found) throw std::invalid_argument("ParamOptimizer: block introuvable: " + name);
        if (ambiguous) throw std::invalid_argument("ParamOptimizer: block '" + name + "' trouvé dans plusieurs scopes.");
        return found;
    }

    void freeze_all_() {
        for (auto& ba : scopes_) if (ba) {
            for (const auto& bn : ba->get_block_names()) ba->at(bn)->freeze();
        }
    }
    void unfreeze_all_() {
        for (auto& ba : scopes_) if (ba) {
            for (const auto& bn : ba->get_block_names()) ba->at(bn)->unfreeze();
        }
    }

    std::vector<std::shared_ptr<BlockAccessor>> scopes_;
    std::vector<Op> ops_;
};
