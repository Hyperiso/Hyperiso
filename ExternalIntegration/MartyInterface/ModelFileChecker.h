#include <string>
#include <fstream>
#include <iostream>
#include <regex>

class ModelFileChecker {
public:
    ModelFileChecker(const std::string& filePath) : filePath(filePath) {}

    bool isAnyModelTemplate() {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filePath);
        }

        std::regex templateRegex(R"(template\s*<[^>]*>)");
        std::regex classRegex(R"(class\s+\w+_Model\s*:\s*public\s+\w+::\w+)");

        std::string line;
        bool foundTemplate = false;

        while (std::getline(file, line)) {
            if (foundTemplate && std::regex_search(line, classRegex)) {
                return true;
            }
            foundTemplate = std::regex_search(line, templateRegex);
        }
        return false;
    }

private:
    std::string filePath;
};