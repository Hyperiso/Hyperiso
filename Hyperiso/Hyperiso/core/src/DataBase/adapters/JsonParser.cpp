#include "JsonParser.h"

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
    LOG_DEBUG("File path:", filename);
    if (!file.is_open()) throw std::runtime_error("Unable to open file for reading");

    std::ostringstream oss;
    oss << file.rdbuf();
    return parse(oss.str());
}