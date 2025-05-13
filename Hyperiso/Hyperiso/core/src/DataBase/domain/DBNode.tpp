#include "DBNode.h"

template <typename... Keys>
Node::Value Node::get(Keys&&... keys) const {
    return getRecursive(data_, std::forward<Keys>(keys)...);
}

template <typename T, typename Key, typename... Rest>
void Node::set(T value, Key&& key, Rest&&... rest) {
    using KeyType = std::decay_t<Key>;
    if constexpr (sizeof...(rest) == 0) {
        if constexpr (std::is_same_v<KeyType, BlockName>)
            data_[std::forward<Key>(key)] = std::forward<T>(value);
        else
            data_[BlockName(std::forward<Key>(key))] = std::forward<T>(value);
    } else {
        auto& node = data_[BlockName(std::forward<Key>(key))];
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
    auto it = map.find(blockKey);
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