#include "DBNode.h"

template <typename... Keys>
Node::Value Node::get(Keys&&... keys) const {
    return getRecursive(data_, std::forward<Keys>(keys)...);
}

template <typename T, typename Key, typename... Rest>
void Node::set(T value, Key&& key, Rest&&... rest) {
    using KeyType = std::decay_t<Key>;
    BlockName blockKey = BlockName(std::forward<Key>(key));

    if constexpr (sizeof...(rest) == 0) {
        auto it = std::find_if(data_.begin(), data_.end(),
            [&](const auto& pair) {
                return pair.first == blockKey;
            });

        if (it != data_.end()) {
            it->second = std::forward<T>(value);
        } else {
            data_[blockKey] = std::forward<T>(value);
        }
    } else {
        auto it = std::find_if(data_.begin(), data_.end(),
            [&](const auto& pair) {
                return pair.first == blockKey;
            });

        if (it == data_.end()) {
            data_[blockKey] = std::make_shared<Node>();
            it = std::find_if(data_.begin(), data_.end(),
                [&](const auto& pair) {
                    return pair.first == blockKey;
                });
        }

        auto& node = it->second;
        if (!std::holds_alternative<std::shared_ptr<Node>>(node)) {
            node = std::make_shared<Node>();
        }
        auto& childNode = std::get<std::shared_ptr<Node>>(node);
        childNode->set(std::forward<T>(value), std::forward<Rest>(rest)...);
    }
}

template <typename Key, typename... Rest>
Node::Value Node::getRecursive(const std::map<BlockName, Value>& map, Key&& key, Rest&&... rest) {
    BlockName blockKey = BlockName(std::forward<Key>(key));
    
    auto it = std::find_if(map.begin(), map.end(),
        [&](const auto& pair) {
            return pair.first == blockKey;
        });
    
    if (it == map.end()) {
        throw std::runtime_error("Key not found");
    }

    if constexpr (sizeof...(rest) == 0) {
        return it->second;
    } else {
        auto node = std::get<std::shared_ptr<Node>>(it->second);
        return node->get(std::forward<Rest>(rest)...);
    }
}