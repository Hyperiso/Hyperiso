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

// template <typename Key, typename... Rest>
// Node::Value Node::getRecursive(const std::map<BlockName, Value>& map, Key&& key, Rest&&... rest) {
//     BlockName blockKey = BlockName(std::forward<Key>(key));
    
//     auto it = std::find_if(map.begin(), map.end(),
//         [&](const auto& pair) {
//             return pair.first == blockKey;
//         });
    
//     if (it == map.end()) {
//         throw std::runtime_error("Key not found");
//     }

//     if constexpr (sizeof...(rest) == 0) {
//         return it->second;
//     } else {
//         auto node = std::get<std::shared_ptr<Node>>(it->second);
//         return node->get(std::forward<Rest>(rest)...);
//     }
// }

static inline bool is_unsigned_number(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    return true;
}

// On suppose que Node::Value est un std::variant<...,
//   std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>>
template <typename Map, typename Key>
Node::Value Node::getRecursive(const Map& map, Key&& key) const {
    BlockName blockKey = BlockName(std::forward<Key>(key));
    auto it = std::find_if(map.begin(), map.end(), [&](const auto& pair){ return pair.first == blockKey; });
    if (it == map.end()) {
        std::ostringstream os;
        os << "Key not found: '" << std::string(blockKey) << "'";
        throw std::runtime_error(os.str());
    }
    return it->second; // FEUILLE => on retourne le variant tel quel
}

template <typename Map, typename Key, typename Next, typename... Rest>
Node::Value Node::getRecursive(const Map& map, Key&& key, Next&& next, Rest&&... rest) const {
    BlockName blockKey = BlockName(std::forward<Key>(key));
    auto it = std::find_if(map.begin(), map.end(), [&](const auto& pair){ return pair.first == blockKey; });
    if (it == map.end()) {
        std::ostringstream os;
        os << "Key not found: '" << std::string(blockKey) << "'";
        throw std::runtime_error(os.str());
    }

    // 1) Cas sous-noeud classique
    if (std::holds_alternative<std::shared_ptr<Node>>(it->second)) {
        auto node = std::get<std::shared_ptr<Node>>(it->second);
        if (!node) {
            std::ostringstream os; os << "Null node at '" << std::string(blockKey) << "'";
            throw std::runtime_error(os.str());
        }
        return node->get(std::forward<Next>(next), std::forward<Rest>(rest)...);
    }

    // 2) Cas liste: on autorise l'accès par index si la "prochaine clé" est un nombre
    if (std::holds_alternative<std::vector<std::shared_ptr<Node>>>(it->second)) {
        const auto& vec = std::get<std::vector<std::shared_ptr<Node>>>(it->second);
        BlockName nextKey = BlockName(std::forward<Next>(next));
        std::string nk = std::string(nextKey);

        if (!is_unsigned_number(nk)) {
            std::ostringstream os;
            os << "Expected numeric index to access list under '" << std::string(blockKey)
               << "', got '" << nk << "'";
            throw std::runtime_error(os.str());
        }
        size_t idx = static_cast<size_t>(std::stoull(nk));
        if (idx >= vec.size()) {
            std::ostringstream os;
            os << "Index " << idx << " out of range for list under '" << std::string(blockKey)
               << "' (size=" << vec.size() << ")";
            throw std::runtime_error(os.str());
        }
        auto node = vec[idx];
        if (!node) {
            std::ostringstream os; os << "Null node at list index " << idx
                                      << " under '" << std::string(blockKey) << "'";
            throw std::runtime_error(os.str());
        }
        if constexpr (sizeof...(Rest) == 0) {
            // On retourne le noeud si c'est la fin du chemin
            return std::static_pointer_cast<Node>(node);
        } else {
            return node->get(std::forward<Rest>(rest)...);
        }
    }

    // 3) Sinon: on ne peut pas descendre au travers d'un scalaire
    {
        std::ostringstream os;
        os << "Cannot descend through non-node value at '" << std::string(blockKey)
           << "' (variant index=" << it->second.index() << ")";
        throw std::runtime_error(os.str());
    }
}