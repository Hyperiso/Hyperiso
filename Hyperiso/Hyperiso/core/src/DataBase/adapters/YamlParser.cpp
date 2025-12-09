#include "YamlParser.h"

bool isInteger(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

bool isDouble(const std::string& s) {
    std::istringstream iss(s);
    double d;
    return (iss >> d) && (iss.eof());
}

bool isBoolean(const std::string& s) {
    return s == "true" || s == "false";
}

bool parseBool(const std::string& s) {
    return s == "true";
}

std::shared_ptr<Node> YAMLParser::parse(const std::string& input) const {
    std::istringstream stream(input);
    return parseYAMLNode(stream, 0);
}

std::shared_ptr<Node> YAMLParser::parseYAMLNode(std::istringstream& stream, int indentLevel) {
    auto node = std::make_shared<Node>();
    std::string line;

    while (std::getline(stream, line)) {
        if(line == ""){
            continue;
        }
        int currentIndent = countLeadingSpaces(line);
        if (currentIndent < indentLevel) {
            stream.seekg(-static_cast<int>(line.size()) - 1, std::ios_base::cur);
            return node;
        }

        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        if (line[0] == '-') {
            processList(node, stream, currentIndent, line);
        } else if (line.find(":") != std::string::npos) {
            processKeyValue(line, node, stream, currentIndent + 2);
        }
    }
    return node;
}

void YAMLParser::processKeyValue(const std::string& line, std::shared_ptr<Node>& node, std::istringstream& stream, int indentLevel) {
    size_t colonPos = line.find(":");
    std::string key = trim(line.substr(0, colonPos));
    std::string value = trim(line.substr(colonPos + 1));
    if (value.empty()) {
        auto childNode = parseYAMLNode(stream, indentLevel);
        node->set(childNode, key);
    } else {
        node->set(parseScalar(value), key);
    }
}

void YAMLParser::processList(std::shared_ptr<Node>& node, std::istringstream& stream, int indentLevel, std::string firstLine) {
    std::map<BlockName, Node::Value> listData;
    size_t index = 0;
    std::string line = firstLine;
    bool first {true};

    Node::Value currentNode;
    auto tempNode = std::make_shared<Node>();

    while (true) {
        if (!first && !std::getline(stream, line)){
            listData[std::to_string(index)] = tempNode->countChildren() != 0 ? tempNode : currentNode;
            break;
        }

        if(trim(line) == ""){
            continue;
        }

        if (!first && countLeadingSpaces(line) < indentLevel) {
            stream.seekg(-static_cast<int>(line.size()) - 1, std::ios_base::cur);
            listData[std::to_string(index)] = tempNode->countChildren() != 0 ? tempNode : currentNode;
            break;
        }

        line = trim(line);
        if (line[0] == '-') {
            if (!first) {
                if (tempNode->countChildren() != 0) {
                    currentNode = tempNode;
                    tempNode = std::make_shared<Node>();
                }
                listData[std::to_string(index)] = currentNode;
                currentNode = Node::Value();
                index++;
            } else {
                first = false;
            }
            line = trim(line.substr(1));
        }

        size_t colonPos = line.find(":");
        
        if (colonPos != std::string::npos) {
            std::string key = trim(line.substr(0, colonPos));
            std::string value = trim(line.substr(colonPos + 1));
            if (value.empty()) {
                auto childNode = parseYAMLNode(stream, indentLevel + 2);
                tempNode->set(childNode, key);
            } else {
                tempNode->set(parseScalar(value), key);
            }
        } else {
            currentNode = parseScalar(trim(line));
        }
    }

    node->setGroup({}, listData);
}



Node::Value YAMLParser::parseScalar(const std::string& value) {
    if (isInteger(value)) {
        return std::stoi(value);
    } else if (isDouble(value)) {
        return std::stod(value);
    } else if (value == "true" || value == "false") {
        return value == "true";
    }
    return value;
}


int YAMLParser::countLeadingSpaces(const std::string& str) {
    return str.find_first_not_of(' ');
}

std::string YAMLParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t");
    return str.substr(first, last - first + 1);
}

bool YAMLParser::isInteger(const std::string& str) {
    return !str.empty() && std::all_of(str.begin(), str.end(), ::isdigit);
}

bool YAMLParser::isDouble(const std::string& str) {
    char* end = nullptr;
    std::strtod(str.c_str(), &end);
    return end != str.c_str() && *end == '\0';
}


bool YAMLParser::isListNode(const std::shared_ptr<Node>& node) const {
    for (const auto& [key, _] : node->getGroup({})) {
        std::string keyStr = key.to_string();

        if (!std::all_of(keyStr.begin(), keyStr.end(), ::isdigit)) {
            return false;
        }
    }
    return true;
}

Node::Value YAMLParser::parseValue(const std::string& value) const {
    if (value.empty()) return "";  

    if (!value.empty() && std::all_of(value.begin(), value.end(), ::isdigit)) {
        return std::stoi(value);
    }

    if (value.find_first_of(".eE") != std::string::npos) {
        try {
            return parseNumber(value);
        } catch (...) {}
    }

    if (value == "true") return true;
    if (value == "false") return false;

    return value;
}

double YAMLParser::parseNumber(const std::string& value) const {
    try {
        return std::stod(value);
    } catch (const std::exception&) {
        throw std::runtime_error("Format de nombre invalide : " + value);
    }
}

void YAMLParser::writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const {
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file for writing: " + filename);
    root->printYAMLToStream(file);
    file << '\n';
    file.flush();
    if (!file)
        throw std::runtime_error("Error while writing YAML file: " + filename);
}

std::shared_ptr<Node> YAMLParser::readFromFile(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

    std::ostringstream oss;
    oss << file.rdbuf();
    return parse(oss.str());
}