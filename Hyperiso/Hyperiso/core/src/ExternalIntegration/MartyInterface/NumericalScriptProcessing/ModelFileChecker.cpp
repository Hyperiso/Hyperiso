#include "ModelFileChecker.h"

#include <algorithm>
#include <cctype>
#include <iterator>
#include <sstream>
#include <unordered_set>

ModelFileChecker::ModelFileChecker(const std::string& filePath) : filePath(filePath) {}

std::string ModelFileChecker::readContents() const {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open MARTY model file: " + filePath);
    }

    std::string contents((std::istreambuf_iterator<char>(file)), {});
    contents.erase(std::remove(contents.begin(), contents.end(), '\r'), contents.end());
    return contents;
}

bool ModelFileChecker::isAnyModelTemplate() const {
    const std::string contents = readContents();
    static const std::regex re(
        R"((?:template\s*<[^>]*>\s*)class\s+[A-Za-z_]\w*\s*(?::|\{|$))",
        std::regex::ECMAScript
    );
    return std::regex_search(contents, re);
}

std::vector<std::string> ModelFileChecker::modelClassCandidates(const std::string& model) {
    auto to_upper = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        return s;
    };
    auto to_lower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return s;
    };

    const std::string upper = to_upper(model);
    const std::string lower = to_lower(model);

    std::vector<std::string> raw = {
        model,
        upper,
        lower,
        model + "_Model",
        upper + "_Model",
        lower + "_Model",
    };

    std::vector<std::string> out;
    std::unordered_set<std::string> seen;
    for (auto& candidate : raw) {
        if (!candidate.empty() && seen.insert(candidate).second) {
            out.emplace_back(std::move(candidate));
        }
    }
    return out;
}

std::string ModelFileChecker::regexEscape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size() * 2);
    for (char c : value) {
        switch (c) {
            case '\\': case '^': case '$': case '.': case '|': case '?':
            case '*': case '+': case '(': case ')': case '[': case ']':
            case '{': case '}':
                escaped.push_back('\\');
                [[fallthrough]];
            default:
                escaped.push_back(c);
        }
    }
    return escaped;
}

bool ModelFileChecker::hasClassDefinition(const std::string& contents,
                                           const std::string& class_name,
                                           bool& is_template) {
    const std::string escaped = regexEscape(class_name);

    const std::regex templated(
        "(?:^|[\\s;{}])template\\s*<[^>]*>\\s*class\\s+" + escaped + R"(\b\s*(?:final\s*)?(?::|\{|$))",
        std::regex::ECMAScript
    );
    if (std::regex_search(contents, templated)) {
        is_template = true;
        return true;
    }

    const std::regex plain(
        "(?:^|[\\s;{}])class\\s+" + escaped + R"(\b\s*(?:final\s*)?(?::|\{|$))",
        std::regex::ECMAScript
    );
    if (std::regex_search(contents, plain)) {
        is_template = false;
        return true;
    }

    return false;
}

ModelClassInfo ModelFileChecker::resolveModelClass(const std::string& model) const {
    const std::string contents = readContents();
    const auto candidates = modelClassCandidates(model);

    for (const auto& candidate : candidates) {
        bool is_template = false;
        if (hasClassDefinition(contents, candidate, is_template)) {
            return ModelClassInfo{candidate, is_template};
        }
    }

    std::ostringstream oss;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (i != 0) {
            oss << ", ";
        }
        oss << candidates[i];
    }

    throw std::runtime_error(
        "Cannot resolve MARTY model class for mty_model_name='" + model +
        "' in file '" + filePath + "'. Tried class names: " + oss.str() + "."
    );
}
