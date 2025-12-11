#include "ModelFileChecker.h"

ModelFileChecker::ModelFileChecker(const std::string& filePath) : filePath(filePath) {}

bool ModelFileChecker::isAnyModelTemplate() const {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filePath);
    }

    std::string contents((std::istreambuf_iterator<char>(file)), {});
    contents.erase(std::remove(contents.begin(), contents.end(), '\r'), contents.end());


    static const std::regex re(
        R"(template\s*<[^>]*>\s*class\s+[A-Za-z_]\w*_Model\s*:\s*public\s+[A-Za-z_]\w*::[A-Za-z_]\w+)",
        std::regex::ECMAScript
    );
    return std::regex_search(contents, re);
}