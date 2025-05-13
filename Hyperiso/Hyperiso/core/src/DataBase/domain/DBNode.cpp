#include "DBNode.h"

Node::Node() = default;

std::vector<BlockName> Node::get_keys() {
    std::vector<BlockName> keys;
    for (auto& [k, _] : this->data_)
        keys.emplace_back(k);
    return keys;
}


std::map<BlockName, Node::Value> Node::getGroup(const std::vector<BlockName>& keys) const {
    const Node* currentNode = this;
    for (const auto& key : keys) {
        auto it = currentNode->data_.find(key);
        if (it == currentNode->data_.end() || !std::holds_alternative<std::shared_ptr<Node>>(it->second)) {
            std::stringstream path;
            for(auto& k : keys) path << k << " ";
            throw std::runtime_error("Key path not found: " + path.str());
        }
        currentNode = std::get<std::shared_ptr<Node>>(it->second).get();
    }
    return currentNode->data_;
}

void Node::setGroup(const std::vector<BlockName>& keys, const std::map<BlockName, Value>& groupData) {
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

        if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
            const auto& list = std::get<std::vector<std::shared_ptr<Node>>>(value);

            bool isFinalList = true;
            for (const auto& node : list) {
                if (!node->contains("")) {
                    isFinalList = false;
                    break;
                }
            }

            if (isFinalList) {
                std::cout << "[ ";
                for (size_t i = 0; i < list.size(); ++i) {
                    list[i]->printValue(list[i]->get(""), level + 4);
                    if (i < list.size() - 1) std::cout << ", ";
                }
                std::cout << " ]";
            }
            else {
                std::cout << "[\n";
                for (size_t i = 0; i < list.size(); ++i) {
                    std::cout << std::string(level + 4, ' ');
                    list[i]->printJSON(level + 4);
                    if (i < list.size() - 1) std::cout << ",";
                    std::cout << "\n";
                }
                std::cout << std::string(level + 2, ' ') << "]";
            }
        } else {
            printValue(value, level);
        }

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
        if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
            std::cout << std::string(level, ' ') << key << ":\n";
            for (const auto& node : std::get<std::vector<std::shared_ptr<Node>>>(value)) {
                std::cout << std::string(level + 2, ' ') << "- ";
                node->printYAML(level + 2);
                std::cout << "\n";
            }
        } else {
            std::cout << std::string(level, ' ') << key << ": ";
            printScalarYAML(value);
            std::cout << "\n";
        }
    }
}


// bool Node::isListNode(const std::shared_ptr<Node>& node) const {
//     for (const auto& [key, _] : node->data_) {
//         if (!std::all_of(key.begin(), key.end(), ::isdigit)) {
//             return false;
//         }
//     }
//     return true;
// }
bool Node::isListNode(const std::shared_ptr<Node>& node) const {
    return std::holds_alternative<std::vector<std::shared_ptr<Node>>>(node->data_.begin()->second);
}




void Node::printValue(const Value& value, int level) const {
    if (std::holds_alternative<BlockName>(value)) {
        std::cout << "\"" << std::get<BlockName>(value) << "\"";
    } else if (std::holds_alternative<int>(value)) {
        std::cout << std::get<int>(value);
    } else if (std::holds_alternative<double>(value)) {
        std::cout << std::get<double>(value);
    } else if (std::holds_alternative<bool>(value)) {
        std::cout << (std::get<bool>(value) ? "true" : "false");
    } else if (std::holds_alternative<std::shared_ptr<Node>>(value)) {
        std::get<std::shared_ptr<Node>>(value)->printJSON(level + 2);
    } else if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
        const auto& list = std::get<std::vector<std::shared_ptr<Node>>>(value);
        std::cout << "[\n";
        for (size_t i = 0; i < list.size(); ++i) {
            std::cout << std::string(level + 2, ' ');
            list[i]->printJSON(level + 2);
            if (i < list.size() - 1) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << std::string(level, ' ') << "]";
    }
}


void Node::printValueToStream(std::ostream& os, const Value& value, int level) const {
    if (std::holds_alternative<BlockName>(value)) {
        os << "\"" << std::get<BlockName>(value) << "\"";
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
    if (std::holds_alternative<BlockName>(value)) {
        std::cout << std::get<BlockName>(value);
    } else if (std::holds_alternative<int>(value)) {
        std::cout << std::get<int>(value);
    } else if (std::holds_alternative<double>(value)) {
        std::cout << std::get<double>(value);
    } else if (std::holds_alternative<bool>(value)) {
        std::cout << (std::get<bool>(value) ? "true" : "false");
    }
}

bool Node::contains(const BlockName& key) const {
    return data_.find(key) != data_.end();
}

int Node::countChildren() const {
    return data_.size();
}