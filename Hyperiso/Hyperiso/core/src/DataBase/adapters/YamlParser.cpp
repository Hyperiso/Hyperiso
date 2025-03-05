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
    auto root = std::make_shared<Node>();

    std::map<int, std::shared_ptr<Node>> indentationMap;
    indentationMap[0] = root;

    std::string line;
    std::string lastKey;
    int lastIndent = 0;
    bool expectingList = false;

    std::vector<std::string> tempStringList;
    std::vector<int> tempIntList;
    std::vector<double> tempDoubleList;
    std::vector<bool> tempBoolList;

    while (std::getline(stream, line)) {
        size_t indent = countLeadingSpaces(line);
        trim(line);

        if (line.empty() || line[0] == '#') continue;

        while (!indentationMap.empty() && indent < lastIndent) {
            indentationMap.erase(lastIndent);
            lastIndent -= 2;
        }

        if (indentationMap.find(indent) == indentationMap.end()) {
            indentationMap[indent] = root;
        }
        auto currentNode = indentationMap[indent];

        if (line.find(":") != std::string::npos) {
            if (expectingList) {
                if (!tempStringList.empty()) currentNode->set(tempStringList, lastKey);
                else if (!tempIntList.empty()) currentNode->set(tempIntList, lastKey);
                else if (!tempDoubleList.empty()) currentNode->set(tempDoubleList, lastKey);
                else if (!tempBoolList.empty()) currentNode->set(tempBoolList, lastKey);

                tempStringList.clear();
                tempIntList.clear();
                tempDoubleList.clear();
                tempBoolList.clear();
                expectingList = false;
            }

            size_t pos = line.find(":");
            std::string key = line.substr(0, pos);
            trim(key);
            std::string value = line.substr(pos + 1);
            trim(value);

            lastKey = key;

            if (value.empty()) {
                expectingList = true;
            } else {
                currentNode->set(parseValue(value), key);
            }

            if (value.empty()) {
                indentationMap[indent + 2] = std::make_shared<Node>();
                currentNode->set(indentationMap[indent + 2], key);
            }
        }
        else if (line.find("- ") == 0) {
            std::string value = line.substr(2);
            trim(value);

            if (lastKey.empty()) {
                throw std::runtime_error("Liste sans clé principale détectée.");
            }

            if (!expectingList) {
                throw std::runtime_error("Liste inattendue sans clé précédente.");
            }

            if (isInteger(value)) tempIntList.push_back(std::stoi(value));
            else if (isDouble(value)) tempDoubleList.push_back(std::stod(value));
            else if (isBoolean(value)) tempBoolList.push_back(parseBool(value));
            else tempStringList.push_back(value);
        }
    }

    if (expectingList && indentationMap.find(lastIndent) != indentationMap.end()) {
        auto lastNode = indentationMap[lastIndent];

        if (!tempStringList.empty()) lastNode->set(tempStringList, lastKey);
        else if (!tempIntList.empty()) lastNode->set(tempIntList, lastKey);
        else if (!tempDoubleList.empty()) lastNode->set(tempDoubleList, lastKey);
        else if (!tempBoolList.empty()) lastNode->set(tempBoolList, lastKey);
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