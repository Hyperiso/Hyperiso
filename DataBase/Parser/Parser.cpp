#include "Parser.h"

std::shared_ptr<Node> JSONParser::parse(const std::string& input) const {
    auto root = std::make_shared<Node>();
    root->set("value", "key1");
    root->set(42, "key2");
    auto child = std::make_shared<Node>();
    child->set(true, "subkey1");
    root->set(child, "child");
    return root;
}

void JSONParser::writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for writing");

    std::ostringstream oss;
    root->printJSONToStream(oss);
    file << oss.str();
    file.close();
}

std::shared_ptr<Node> JSONParser::readFromFile(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

    std::ostringstream oss;
    oss << file.rdbuf();
    return parse(oss.str());
}

std::shared_ptr<Node> YAMLParser::parse(const std::string& input) const {
    auto root = std::make_shared<Node>();
    root->set("yamlKey", "yamlValue");
    return root;
}

void YAMLParser::writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for writing");
    root->printYAML();
    file.close();
}

std::shared_ptr<Node> YAMLParser::readFromFile(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

    std::ostringstream oss;
    oss << file.rdbuf();
    return parse(oss.str());
}

std::unique_ptr<Parser> ParserFactory::createParser(Type type) {
    switch (type) {
    case Type::JSON:
        return std::make_unique<JSONParser>();
    case Type::YAML:
        return std::make_unique<YAMLParser>();
    default:
        throw std::invalid_argument("Unknown parser type");
    }
}