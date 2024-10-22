#include <string>
#include <fstream>
#include <iostream>
#include <regex>

class ModelFileChecker {
public:
    ModelFileChecker(const std::string& filePath) : filePath(filePath) {}

    bool isTemplateClass(const std::string& className) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filePath);
        }

        std::string pattern = "template\\s*<[^>]*>\\s*class\\s+" + className;
        std::regex templateClassRegex(pattern);

        std::string line;
        while (std::getline(file, line)) {
            if (std::regex_search(line, templateClassRegex)) {
                return true;
            }
        }
        return false;
    }

    bool isAnyModelTemplate() {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filePath);
        }

        std::regex templateClassRegex(R"(template\s*<[^>]*>\s*class\s+\w*_Model)");

        std::string line;
        while (std::getline(file, line)) {
            if (std::regex_search(line, templateClassRegex)) {
                return true;
            }
        }
        return false;
    }

private:
    std::string filePath;
};