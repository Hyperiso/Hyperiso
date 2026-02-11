#include "YamlInputReader.h"

YamlInputReader::YamlInputReader(const std::string& filename) {
    root_ = yaml_p.readFromFile(filename);
    if (!root_) {
        throw std::runtime_error("YamlInputReader: root DBNode is null after parsing file: " + filename);
    }
}

std::string YamlInputReader::unquote(std::string s) {
    if (s.size() >= 2) {
        char a = s.front();
        char b = s.back();
        if ((a == '"' && b == '"') || (a == '\'' && b == '\'')) {
            return s.substr(1, s.size() - 2);
        }
    }
    return s;
}

std::string YamlInputReader::as_string(const DBNode::Value& v, const std::string& ctx) {
    if (std::holds_alternative<BlockName>(v)) {
        return std::get<BlockName>(v).to_string();
    }

    throw std::runtime_error("YamlInputReader: expected string(BlockName) for " + ctx);
}

double YamlInputReader::as_double(const DBNode::Value& v, const std::string& ctx) {
    if (std::holds_alternative<double>(v)) return std::get<double>(v);
    if (std::holds_alternative<int>(v))    return static_cast<double>(std::get<int>(v));
    throw std::runtime_error("YamlInputReader: expected number(int/double) for " + ctx);
}

std::shared_ptr<DBNode> YamlInputReader::as_node(const DBNode::Value& v, const std::string& ctx) {
    if (std::holds_alternative<std::shared_ptr<DBNode>>(v)) {
        return std::get<std::shared_ptr<DBNode>>(v);
    }
    throw std::runtime_error("YamlInputReader: expected DBNode for " + ctx);
}

std::vector<YamlInputParam> YamlInputReader::get_input_params() const {
    std::vector<YamlInputParam> out;

    auto pSpecsNode = as_node(root_->get("p_specs"), "p_specs");
    auto items = pSpecsNode->getGroup({}); 

    out.reserve(items.size());

    for (const auto& [idxKey, itemVal] : items) {
        (void)idxKey;

        auto itemNode = as_node(itemVal, "p_specs item");

        auto block_raw = as_string(itemNode->get("block"), "p_specs.block");
        auto code_raw  = as_string(itemNode->get("code"),  "p_specs.code");

        std::string block = unquote(block_raw);
        std::string code  = unquote(code_raw);

        out.push_back(YamlInputParam{
            .block_name = block,
            .pdg_code   = LhaID(code)
        });
    }

    return out;
}

std::vector<YamlScanParam> YamlInputReader::get_scan_params() const {
    std::vector<YamlScanParam> out;

    auto scanNode = as_node(root_->get("scan_params"), "scan_params");
    auto items = scanNode->getGroup({});

    out.reserve(items.size());

    for (const auto& [idxKey, itemVal] : items) {
        (void)idxKey;

        auto itemNode = as_node(itemVal, "scan_params item");

        auto block_raw = as_string(itemNode->get("block"), "scan_params.block");
        auto code_raw  = as_string(itemNode->get("code"),  "scan_params.code");

        double minv = as_double(itemNode->get("min_val"), "scan_params.min_val");
        double maxv = as_double(itemNode->get("max_val"), "scan_params.max_val");
        double step = as_double(itemNode->get("step"),    "scan_params.step");

        std::string block = unquote(block_raw);
        std::string code  = unquote(code_raw);

        out.push_back(YamlScanParam{
            .block_name = block,
            .pdg_code   = LhaID(code),
            .min_val    = minv,
            .max_val    = maxv,
            .step_val   = step
        });
    }

    return out;
}
