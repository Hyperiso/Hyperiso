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
        auto it = std::find_if(currentNode->data_.begin(), currentNode->data_.end(),
        [&](const auto& pair) {
            return pair.first == key;
        });
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
        auto it = std::find_if(currentNode->data_.begin(), currentNode->data_.end(),
        [&](const auto& pair) {
            return pair.first == key;
        });

    if (it == currentNode->data_.end()) {
        currentNode->data_[key] = std::make_shared<Node>();
        it = std::find_if(currentNode->data_.begin(), currentNode->data_.end(),
            [&](const auto& pair) {
                return pair.first == key;
            });
    }

    auto& value = it->second;
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
        const auto& key = it->first;
        const auto& value = it->second;

        os << std::string(level + 2, ' ') << "\"" << key << "\": ";

        if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
            const auto& list = std::get<std::vector<std::shared_ptr<Node>>>(value);

            bool isFinalList = true;
            for (const auto& node : list) {
                if (!node || !node->contains("")) { isFinalList = false; break; }
            }

            if (isFinalList) {
                os << "[ ";
                for (size_t i = 0; i < list.size(); ++i) {
                    printValueToStream(os, list[i]->get(""), level + 4);
                    if (i + 1 < list.size()) os << ", ";
                }
                os << " ]";
            } else {
                os << "[\n";
                for (size_t i = 0; i < list.size(); ++i) {
                    os << std::string(level + 4, ' ');
                    list[i]->printJSONToStream(os, level + 4);
                    if (i + 1 < list.size()) os << ",";
                    os << "\n";
                }
                os << std::string(level + 2, ' ') << "]";
            }
        } else {
            printValueToStream(os, value, level);
        }

        if (std::next(it) != data_.end()) os << ",";
        os << "\n";
    }
    os << std::string(level, ' ') << "}";
}

void Node::printYAML(int level) const {
    for (const auto& [key, value] : data_) {
        if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
            const auto& list = std::get<std::vector<std::shared_ptr<Node>>>(value);
            std::cout << std::string(level, ' ') << key << ":\n";

            bool isFinalList = true;
            for (const auto& node : list) {
                if (!node || !node->contains("")) { isFinalList = false; break; }
            }

            if (isFinalList) {
                for (const auto& node : list) {
                    std::cout << std::string(level + 2, ' ') << "- ";
                    printScalarYAML(node->get(""));
                    std::cout << "\n";
                }
            } else {
                for (const auto& node : list) {
                    std::cout << std::string(level + 2, ' ') << "- ";
                    if (node) {
                        std::cout << "\n";
                        node->printYAML(level + 4);
                    } else {
                        std::cout << "null\n";
                    }
                }
            }
            continue;
        }

        if (std::holds_alternative<std::shared_ptr<Node>>(value)) {
            auto child = std::get<std::shared_ptr<Node>>(value);
            std::cout << std::string(level, ' ') << key << ":\n";
            if (child) child->printYAML(level + 2);
            else       std::cout << std::string(level + 2, ' ') << "null\n";
            continue;
        }

        std::cout << std::string(level, ' ') << key << ": ";
        printScalarYAML(value);
        std::cout << "\n";
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
        auto ptr = std::get<std::shared_ptr<Node>>(value);
        if (ptr) ptr->printJSONToStream(os, level + 2);
        else     os << "null";
    } else if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
        const auto& list = std::get<std::vector<std::shared_ptr<Node>>>(value);
        os << "[\n";
        for (size_t i = 0; i < list.size(); ++i) {
            os << std::string(level + 2, ' ');
            if (list[i]) list[i]->printJSONToStream(os, level + 2);
            else         os << "null";
            if (i + 1 < list.size()) os << ",";
            os << "\n";
        }
        os << std::string(level, ' ') << "]";
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
    return std::any_of(data_.begin(), data_.end(),
        [&](const auto& pair) {
            return pair.first == key;
        });
}

int Node::countChildren() const {
    return data_.size();
}

void Node::printScalarYAMLToStream(std::ostream& os, const Value& value) const {
    if (std::holds_alternative<BlockName>(value)) {
        os << std::get<BlockName>(value);
    } else if (std::holds_alternative<int>(value)) {
        os << std::get<int>(value);
    } else if (std::holds_alternative<double>(value)) {
        os << std::get<double>(value);
    } else if (std::holds_alternative<bool>(value)) {
        os << (std::get<bool>(value) ? "true" : "false");
    }
}

void Node::printYAMLToStream(std::ostream& os, int level) const {
    for (const auto& [key, value] : data_) {
        if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(value)) {
            const auto& list = std::get<std::vector<std::shared_ptr<Node>>>(value);
            os << std::string(level, ' ') << key << ":\n";

            bool isFinalList = true;
            for (const auto& node : list) {
                if (!node || !node->contains("")) { isFinalList = false; break; }
            }

            if (isFinalList) {
                for (const auto& node : list) {
                    os << std::string(level + 2, ' ') << "- ";
                    printScalarYAMLToStream(os, node->get(""));
                    os << "\n";
                }
            } else {
                for (const auto& node : list) {
                    os << std::string(level + 2, ' ') << "- ";
                    if (node) {
                        os << "\n";
                        node->printYAMLToStream(os, level + 4);
                    } else {
                        os << "null\n";
                    }
                }
            }
            continue;
        }

        if (std::holds_alternative<std::shared_ptr<Node>>(value)) {
            auto child = std::get<std::shared_ptr<Node>>(value);
            os << std::string(level, ' ') << key << ":\n";
            if (child) child->printYAMLToStream(os, level + 2);
            else       os << std::string(level + 2, ' ') << "null\n";
            continue;
        }

        os << std::string(level, ' ') << key << ": ";
        printScalarYAMLToStream(os, value);
        os << "\n";
    }
}