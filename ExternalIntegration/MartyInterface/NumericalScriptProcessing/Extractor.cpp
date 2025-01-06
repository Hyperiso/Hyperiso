#include "Extractor.h"
#include <regex>
#include <fstream>

std::vector<Extractor::Parameter> Extractor::extract(const std::string& filename) {
    std::vector<Parameter> params;
    std::ifstream file(filename);
    std::string line;
    std::regex paramRegex(R"(csl::InitSanitizer<real_t>\s+(\w+)\s+\{\s+\"(\w+)\"\s+\};)");

    while (std::getline(file, line)) {
        std::smatch match;
        if (std::regex_search(line, match, paramRegex)) {
            Parameter param{match[1], match[2], false};
            params.push_back(param);
        }
    }

    std::ifstream file_2(filename);
    std::string line_2;
    std::regex paramRegex_2(R"(csl::InitSanitizer<complex_t>\s+(\w+)\s+\{\s+\"(\w+)\"\s+\};)");

    while (std::getline(file_2, line_2)) {
        std::smatch match;
        if (std::regex_search(line_2, match, paramRegex_2)) {
            Parameter param{match[1], match[2], true};
            params.push_back(param);
        }
    }

    return params;
}
