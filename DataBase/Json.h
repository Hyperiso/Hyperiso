#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <variant>
#include <initializer_list>
#include <stdexcept>

class Node {
public:
    using Value = std::variant<std::string, int, double, bool, std::shared_ptr<Node>>;

private:
    std::map<std::string, Value> data_;

public:
    Node() = default;

    template <typename... Keys>
    Value get(Keys&&... keys) const {
        return getRecursive(data_, std::forward<Keys>(keys)...);
    }

    template <typename T, typename Key, typename... Rest>
    void set(T value, Key&& key, Rest&&... rest) {
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

    void printJSON(int level = 0) const {
        std::cout << "{\n";
        for (const auto& [key, value] : data_) {
            std::cout << std::string(level + 2, ' ') << "\"" << key << "\": ";
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
            std::cout << ",\n";
        }
        std::cout << std::string(level, ' ') << "}";
    }

    void printJSONToStream(std::ostream& os, int level = 0) const {
        os << "{\n";
        for (const auto& [key, value] : data_) {
            os << std::string(level + 2, ' ') << "\"" << key << "\": ";
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
            os << ",\n";
        }
        os << std::string(level, ' ') << "}";
    }
    

private:
    template <typename Key, typename... Rest>
    static Value getRecursive(const std::map<std::string, Value>& map, Key&& key, Rest&&... rest) {
        auto it = map.find(std::forward<Key>(key));
        if (it == map.end())
            throw std::runtime_error("Key not found");
        if constexpr (sizeof...(rest) == 0) {
            return it->second;
        } else {
            auto node = std::get<std::shared_ptr<Node>>(it->second);
            return node->get(std::forward<Rest>(rest)...);
        }
    }
};



class Parser {
public:
    virtual ~Parser() = default;

    virtual std::shared_ptr<Node> parse(const std::string& input) const = 0;
    virtual void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const = 0;
    virtual std::shared_ptr<Node> readFromFile(const std::string& filename) const = 0;
};

class JSONParser : public Parser {
public:
    std::shared_ptr<Node> parse(const std::string& input) const override {
        auto root = std::make_shared<Node>();
        root->set("value", "key1");
        root->set(42, "key2");
        auto child = std::make_shared<Node>();
        child->set(true, "subkey1");
        root->set(child, "child");
        return root;
    }

    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override {
        std::ofstream file(filename);
        if (!file.is_open()) throw std::runtime_error("Unable to open file for writing");

        std::ostringstream oss;
        root->printJSONToStream(oss);
        file << oss.str();
        file.close();
    }


    std::shared_ptr<Node> readFromFile(const std::string& filename) const override {
        std::ifstream file(filename);
        if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

        std::ostringstream oss;
        oss << file.rdbuf();
        return parse(oss.str());
    }
};

class YAMLParser : public Parser {
public:
    std::shared_ptr<Node> parse(const std::string& input) const override {
        auto root = std::make_shared<Node>();
        root->set("yamlKey", "yamlValue");
        return root;
    }

    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override {
        std::ofstream file(filename);
        if (!file.is_open()) throw std::runtime_error("Unable to open file for writing");
        root->printJSON();
        file.close();
    }

    std::shared_ptr<Node> readFromFile(const std::string& filename) const override {
        std::ifstream file(filename);
        if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

        std::ostringstream oss;
        oss << file.rdbuf();
        return parse(oss.str());
    }
};

class ParserFactory {
public:
    enum class Type { JSON, YAML };

    static std::unique_ptr<Parser> createParser(Type type) {
        switch (type) {
        case Type::JSON:
            return std::make_unique<JSONParser>();
        case Type::YAML:
            return std::make_unique<YAMLParser>();
        default:
            throw std::invalid_argument("Unknown parser type");
        }
    }
};
