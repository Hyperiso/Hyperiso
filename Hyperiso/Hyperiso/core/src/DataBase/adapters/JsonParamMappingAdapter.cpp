#include "JsonParamMappingAdapter.h"
#include "Logger.h"

std::unordered_map<std::string, InterpretedParam>
JsonParamMappingAdapter::loadFromFile(const std::string& filePath) const {
    std::unordered_map<std::string, InterpretedParam> out;

    if (!parser_) {
        throw std::runtime_error("JsonParamMappingAdapter: parser nul");
    }

    auto root = parser_->readFromFile(filePath);
    if (!root) {
        throw std::runtime_error("JsonParamMappingAdapter: Node root nul");
    }

    auto keys = root->get_keys();
    for (const auto& paramName : keys) {
        try {
            const auto blockV = root->get(paramName, "block");
            const auto codeV  = root->get(paramName, "pdgCode");

            const std::string block = asString(blockV);
            const std::string code  = asString(codeV);

            LhaID lh = LhaID(code);

            out[static_cast<std::string>(paramName)] = InterpretedParam{ block, lh };
        } catch (const std::exception& e) {
            LOG_ERROR("JSON mapping", std::string("Key '") + static_cast<std::string>(paramName) + "': " + e.what());
        }
    }
    return out;
}
