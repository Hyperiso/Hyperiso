#include "JsonParameters.h"

std::string JSONBlock::toJSON(int indentLevel) const {
    std::stringstream ss;
    std::string indent(indentLevel * 2, ' ');
    ss << "{\n";
    bool first = true;
    for (const auto& [pdgCode, value] : pdgCodeValues) {
        if (!first) ss << ",\n";
        ss << indent << "\"" << pdgCode << "\": " << value;
        first = false;
    }
    ss << "\n" << std::string((indentLevel - 1) * 2, ' ') << "}";
    return ss.str();
}

void JSONBlock::fromJSON(const std::string& json) {
    pdgCodeValues.clear();
    std::istringstream stream(json);
    char c;
    int pdgCode;
    double value;
    while (stream >> c) {
        if (c == '"') {
            stream >> pdgCode;
            stream >> c >> c;
            stream >> value;
            pdgCodeValues[pdgCode] = value;
            stream >> c;
            if (c == '}') break;
        }
    }
}

std::string JSONParser::toJSON(int indentLevel) const {
    std::stringstream ss;
    std::string indent(indentLevel * 2, ' ');
    ss << "{\n";
    bool first = true;
    for (const auto& [blockName, block] : blocks) {
        if (!first) ss << ",\n";
        ss << indent << "\"" << blockName << "\": " << block.toJSON(indentLevel + 1);
        first = false;
    }
    ss << "\n" << std::string(indentLevel * 2, ' ') << "}";
    return ss.str();
}

void JSONParser::fromJSON(const std::string& json) {
    blocks.clear();
    std::istringstream stream(json);
    char c;
    stream >> c;
    while (stream >> c && c != '}') {
        if (c == '"') {
            std::string blockName;
            std::getline(stream, blockName, '"');
            stream >> c >> c;
            std::string blockData;
            std::getline(stream, blockData, '}');
            blockData += '}';
            JSONBlock block;
            block.fromJSON(blockData);
            blocks[blockName] = block;
            stream >> c;
        }
    }
}

std::unordered_map<int, JSONParser*> JSONParser::instances;