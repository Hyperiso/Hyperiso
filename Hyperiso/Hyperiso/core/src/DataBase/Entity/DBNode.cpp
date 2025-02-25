#include "DBNode.h"

Node::Node() = default;

std::vector<std::string> Node::get_keys() {
    std::vector<std::string> keys;
    for (auto& [k, _] : this->data_)
        keys.emplace_back(k);
    return keys;
}

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

std::map<std::string, Node::Value> Node::getGroup(const std::vector<std::string>& keys) const {
    const Node* currentNode = this;
    for (const auto& key : keys) {
        auto it = currentNode->data_.find(key);
        if (it == currentNode->data_.end() || !std::holds_alternative<std::shared_ptr<Node>>(it->second)) {
            throw std::runtime_error("Key path not found");
        }
        currentNode = std::get<std::shared_ptr<Node>>(it->second).get();
    }
    return currentNode->data_;
}

void Node::setGroup(const std::vector<std::string>& keys, const std::map<std::string, Value>& groupData) {
    Node* currentNode = this;
    for (const auto& key : keys) {
        auto& value = currentNode->data_[key];
        if (!std::holds_alternative<std::shared_ptr<Node>>(value)) {
            value = std::make_shared<Node>();
        }
        currentNode = std::get<std::shared_ptr<Node>>(value).get();
    }
    currentNode->data_ = groupData;
}

void Node::printJSON(int level) const {
    std::cout << "{\n";
    for (auto it = data_.begin(); it != data_.end(); ++it) {
        const auto& [key, value] = *it;
        std::cout << std::string(level + 2, ' ') << "\"" << key << "\": ";
        printValue(value, level);
        if (std::next(it) != data_.end()) {
            std::cout << ",";
        }
        std::cout << "\n";
    }
    std::cout << std::string(level, ' ') << "}";
}

void Node::printJSONToStream(std::ostream& os, int level) const {
    os << "{\n";
    for (auto it = data_.begin(); it != data_.end(); ++it) {
        const auto& [key, value] = *it;
        os << std::string(level + 2, ' ') << "\"" << key << "\": ";
        printValueToStream(os, value, level);
        if (std::next(it) != data_.end()) {
            os << ",";
        }
        os << "\n";
    }
    os << std::string(level, ' ') << "}";
}

void Node::printYAML(int level) const {
    for (const auto& [key, value] : data_) {
        std::cout << std::string(level, ' ') << key << ": ";
        if (std::holds_alternative<std::string>(value) || std::holds_alternative<int>(value) ||
            std::holds_alternative<double>(value) || std::holds_alternative<bool>(value)) {
            printScalarYAML(value);
            std::cout << "\n";
        } else if (std::holds_alternative<std::shared_ptr<Node>>(value)) {
            std::cout << "\n";
            std::get<std::shared_ptr<Node>>(value)->printYAML(level + 2);
        }
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

void Node::printValue(const Value& value, int level) const {
    if (std::holds_alternative<std::string>(value)) {
        std::cout << "\"" << std::get<std::string>(value) << "\"";
    } else if (std::holds_alternative<int>(value)) {
        std::cout << std::get<int>(value);
    } else if (std::holds_alternative<double>(value)) {
        std::cout << std::get<double>(value);
    } else if (std::holds_alternative<bool>(value)) {
        std::cout << (std::get<bool>(value) ? "true" : "false");
    } else if (std::holds_alternative<std::shared_ptr<Node>>(value)) {
        std::get<std::shared_ptr<Node>>(value)->printJSON(level + 2);
    }
}

void Node::printValueToStream(std::ostream& os, const Value& value, int level) const {
    if (std::holds_alternative<std::string>(value)) {
        os << "\"" << std::get<std::string>(value) << "\"";
    } else if (std::holds_alternative<int>(value)) {
        os << std::get<int>(value);
    } else if (std::holds_alternative<double>(value)) {
        os << std::get<double>(value);
    } else if (std::holds_alternative<bool>(value)) {
        os << (std::get<bool>(value) ? "true" : "false");
    } else if (std::holds_alternative<std::shared_ptr<Node>>(value)) {
        std::get<std::shared_ptr<Node>>(value)->printJSONToStream(os, level + 2);
    }
}

void Node::printScalarYAML(const Value& value) const {
    if (std::holds_alternative<std::string>(value)) {
        std::cout << std::get<std::string>(value);
    } else if (std::holds_alternative<int>(value)) {
        std::cout << std::get<int>(value);
    } else if (std::holds_alternative<double>(value)) {
        std::cout << std::get<double>(value);
    } else if (std::holds_alternative<bool>(value)) {
        std::cout << (std::get<bool>(value) ? "true" : "false");
    }
}