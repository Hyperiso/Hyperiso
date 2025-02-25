#include "Parser.h"

std::shared_ptr<Node> JSONParser::parse(const std::string& input) const {
    size_t index = 0;
    return parseObject(input, index);

}

std::shared_ptr<Node> JSONParser::parseObject(const std::string& input, size_t& index) const {
    auto root = std::make_shared<Node>();
    if (input[index] != '{') throw std::runtime_error("Invalid JSON");
    ++index;

    while (index < input.size()) {
        skipWhitespace(input, index);
        if (input[index] == '}') { ++index; break; }
        
        std::string key = parseString(input, index);
        skipWhitespace(input, index);
        if (input[index] != ':') throw std::runtime_error("Expected ':' in JSON");
        ++index;
        skipWhitespace(input, index);
        
        root->set(parseValue(input, index), key);
        skipWhitespace(input, index);
        if (input[index] == ',') ++index;
    }
    return root;
}

std::string JSONParser::parseString(const std::string& input, size_t& index) const {
    if (input[index] != '\"') throw std::runtime_error("Expected '\"' in JSON");
    ++index;
    std::string result;
    while (index < input.size() && input[index] != '\"') {
        result += input[index++];
    }
    if (input[index] != '\"') throw std::runtime_error("Expected closing '\"'");
    ++index;
    return result;
}

Node::Value JSONParser::parseValue(const std::string& input, size_t& index) const {
    skipWhitespace(input, index);

    if (input[index] == '{') return parseObject(input, index);
    if (input[index] == '[') return parseArray(input, index);
    if (input[index] == '\"') return parseString(input, index);
    if (isdigit(input[index]) || input[index] == '-') return parseNumber(input, index);
    if (input.substr(index, 4) == "true") { index += 4; return true; }
    if (input.substr(index, 5) == "false") { index += 5; return false; }

    throw std::runtime_error("Invalid value in JSON");
}


std::shared_ptr<Node> JSONParser::parseArray(const std::string& input, size_t& index) const {
    auto listNode = std::make_shared<Node>();
    ++index; // Skip '['

    int elementIndex = 0;
    while (index < input.size()) {
        skipWhitespace(input, index);
        if (input[index] == ']') { ++index; break; }

        listNode->set(parseValue(input, index), std::to_string(elementIndex++));
        skipWhitespace(input, index);

        if (input[index] == ',') ++index;
    }
    return listNode;
}

double JSONParser::parseNumber(const std::string& input, size_t& index) const {
    std::string num;
    bool hasDecimal = false;
    bool hasExponent = false;

    while (index < input.size() && (isdigit(input[index]) || input[index] == '.' || input[index] == 'e' || input[index] == 'E' || input[index] == '-' || input[index] == '+')) {
        if (input[index] == '.') {
            if (hasDecimal || hasExponent) break;
            hasDecimal = true;
        } else if (input[index] == 'e' || input[index] == 'E') {
            if (hasExponent) break;
            hasExponent = true;
            if (index + 1 < input.size() && (input[index + 1] == '-' || input[index + 1] == '+')) {
                num += input[index++];
            }
        }
        num += input[index++];
    }
    try {
        return std::stod(num);
    } catch (...) {
        throw std::runtime_error("Invalid number format in JSON");
    }
}

void JSONParser::skipWhitespace(const std::string& input, size_t& index) const {
    while (index < input.size() && isspace(input[index])) ++index;
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
    LOG_INFO("File path:", filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

    std::ostringstream oss;
    oss << file.rdbuf();
    return parse(oss.str());
}

std::shared_ptr<Node> YAMLParser::parse(const std::string& input) const {
    std::istringstream stream(input);
    auto root = std::make_shared<Node>();

    std::stack<std::shared_ptr<Node>> nodeStack;
    std::map<int, std::shared_ptr<Node>> indentationMap;
    nodeStack.push(root);
    indentationMap[0] = root;

    std::string line;
    std::string lastKey;
    int lastIndent = 0;
    bool lastWasList = false;

    while (std::getline(stream, line)) {
        size_t indent = countLeadingSpaces(line);
        trim(line);

        if (line.empty() || line[0] == '#') continue;

        if (indent < lastIndent) {
            while (!indentationMap.empty() && indent < lastIndent) {
                indentationMap.erase(lastIndent);
                lastIndent -= 2;
            }
        }
        lastIndent = indent;

        auto currentNode = indentationMap[indent];

        if (line.find(":") != std::string::npos) {
            size_t pos = line.find(":");
            std::string key = line.substr(0, pos);
            trim(key);
            std::string value = line.substr(pos + 1);
            trim(value);

            lastKey = key;
            lastWasList = false;

            if (value.empty()) {
                if (!currentNode->contains(key)) {
                    currentNode->set(std::make_shared<Node>(), key);
                }
                auto newNode = std::get<std::shared_ptr<Node>>(currentNode->get(key));
                nodeStack.push(newNode);
                indentationMap[indent + 2] = newNode;
            } else {
                currentNode->set(parseValue(value), key);
            }
        } else if (line.find("- ") == 0) {
            std::string value = line.substr(2);
            trim(value);

            if (lastKey.empty()) {
                throw std::runtime_error("Liste sans clé principale détectée.");
            }

            if (!currentNode->contains(lastKey) || !isListNode(std::get<std::shared_ptr<Node>>(currentNode->get(lastKey)))) {
                currentNode->set(std::make_shared<Node>(), lastKey);
            }

            auto listNode = std::get<std::shared_ptr<Node>>(currentNode->get(lastKey));

            int index = listNode->countChildren();
            listNode->set(parseValue(value), std::to_string(index));

            lastWasList = true;
        }
    }
    return root;
}

bool YAMLParser::isListNode(const std::shared_ptr<Node>& node) const {
    for (const auto& [key, _] : node->getGroup({})) {
        if (!std::all_of(key.begin(), key.end(), ::isdigit)) {
            return false;
        }
    }
    return true;
}

size_t YAMLParser::countLeadingSpaces(const std::string& line) const {
    size_t count = 0;
    while (count < line.size() && isspace(line[count])) count++;
    return count;
}

void YAMLParser::trim(std::string& str) const {
    while (!str.empty() && isspace(str.front())) str.erase(str.begin());
    while (!str.empty() && isspace(str.back())) str.pop_back();
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