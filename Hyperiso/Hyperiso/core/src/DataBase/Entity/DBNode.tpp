#include "DBNode.h"

template <typename... Keys>
Node::Value Node::get(Keys&&... keys) const {
    return getRecursive(data_, std::forward<Keys>(keys)...);
}

template <typename T, typename Key, typename... Rest>
void Node::set(T value, Key&& key, Rest&&... rest) {
    if constexpr (sizeof...(rest) == 0) {
        data_[std::string(std::forward<Key>(key))] = std::forward<T>(value);
    } else {
        auto& node = data_[std::string(std::forward<Key>(key))];
        if (!std::holds_alternative<std::shared_ptr<Node>>(node)) {
            node = std::make_shared<Node>();
        }
        auto& childNode = std::get<std::shared_ptr<Node>>(node);
        childNode->set(std::forward<T>(value), std::forward<Rest>(rest)...);
    }
}

template <typename Key, typename... Rest>
Node::Value Node::getRecursive(const std::map<std::string, Value>& map, Key&& key, Rest&&... rest) {
    auto it = map.find(std::forward<Key>(key));
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