#include "MappingDataBase.h"
#include <fstream>
#include <iostream>
#include <sstream>

std::unordered_map<std::string, std::shared_ptr<MappingDatabase>> MappingDatabase::instances;

std::shared_ptr<MappingDatabase> MappingDatabase::getInstance(const std::string& instanceName, const std::string& jsonFilePath) {
    if (instances.find(instanceName) != instances.end()) {
        return instances[instanceName];
    }
    
    if (!jsonFilePath.empty()) {
        instances[instanceName] = std::shared_ptr<MappingDatabase>(new MappingDatabase(jsonFilePath));
    } else {
        std::cerr << "Instance " << instanceName << " not found and no file path provided to create it." << std::endl;
    }
    return instances[instanceName];
}

MappingDatabase::MappingDatabase(const std::string& jsonFilePath) {
    loadFromJson(jsonFilePath);
}

void MappingDatabase::loadFromJson(const std::string& jsonFilePath) {
    std::ifstream file(jsonFilePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open JSON file: " << jsonFilePath << std::endl;
        return;
    }

    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    size_t pos = 0;
    while ((pos = content.find("\"", pos)) != std::string::npos) {
        size_t start = pos + 1;
        size_t end = content.find("\"", start);
        std::string key = content.substr(start, end - start);
        pos = end + 1;

        size_t blockPos = content.find("\"block\"", pos);
        if (blockPos == std::string::npos) break;
        start = content.find("\"", blockPos + 8) + 1;
        end = content.find("\"", start);
        std::string block = content.substr(start, end - start);
        pos = end + 1;

        size_t pdgPos = content.find("\"pdgCode\"", pos);
        if (pdgPos == std::string::npos) break;
        start = content.find(":", pdgPos) + 1;
        end = content.find_first_of(",}", start);
        std::string codeStr = content.substr(start, end - start);

        codeStr.erase(0, codeStr.find_first_not_of(" \t"));
        codeStr.erase(codeStr.find_last_not_of(" ,\t") + 1);
        int pdgCode = 0;
        try {
            pdgCode = std::stoi(codeStr);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid integer format for pdgCode: " << codeStr << std::endl;
            pdgCode = -1;
        } catch (const std::out_of_range& e) {
            std::cerr << "Integer out of range for pdgCode: " << codeStr << std::endl;
            pdgCode = -1;
        }

        if (!key.empty() && !block.empty()) {
            paramsMap[key] = {block, pdgCode};
        }

        pos = content.find("}", end) + 1;
    }
}

const std::unordered_map<std::string, InterpretedParam>& MappingDatabase::getParams() const {
    return paramsMap;
}
