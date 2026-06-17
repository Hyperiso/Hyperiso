#include "NuisanceReader.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <sstream>

namespace {

std::string strip_optional_quotes(std::string s) {
    if (s.size() >= 2) {
        const char first = s.front();
        const char last  = s.back();

        if ((first == '"'  && last == '"') ||
            (first == '\'' && last == '\'')) {
            return s.substr(1, s.size() - 2);
        }
    }
    return s;
}

bool is_numeric_string(const std::string& s) {
    if (s.empty()) return false;
    char* end = nullptr;
    std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

bool looks_like_nuisance_entry(const std::shared_ptr<DBNode>& node) {
    return node
        && node->contains("block")
        && node->contains("code");
}

const char* value_kind(const DBNode::Value& value) {
    if (std::holds_alternative<BlockName>(value)) return "string";
    if (std::holds_alternative<int>(value)) return "int";
    if (std::holds_alternative<double>(value)) return "double";
    if (std::holds_alternative<bool>(value)) return "bool";
    if (std::holds_alternative<std::shared_ptr<DBNode>>(value)) return "node";
    if (std::holds_alternative<std::vector<std::shared_ptr<DBNode>>>(value)) return "list";
    return "unknown";
}

void collect_entry_nodes_from_value(const DBNode::Value& value,
                                    std::vector<std::shared_ptr<DBNode>>& out);

void collect_entry_nodes_from_node(const std::shared_ptr<DBNode>& node,
                                   std::vector<std::shared_ptr<DBNode>>& out)
{
    if (!node) {
        return;
    }

    if (looks_like_nuisance_entry(node)) {
        out.push_back(node);
        return;
    }

    for (const auto& key : node->get_keys()) {
        const auto child = node->get(key);
        collect_entry_nodes_from_value(child, out);
    }
}

void collect_entry_nodes_from_value(const DBNode::Value& value,
                                    std::vector<std::shared_ptr<DBNode>>& out)
{
    if (std::holds_alternative<std::vector<std::shared_ptr<DBNode>>>(value)) {
        const auto& entries = std::get<std::vector<std::shared_ptr<DBNode>>>(value);
        for (const auto& entry : entries) {
            collect_entry_nodes_from_node(entry, out);
        }
        return;
    }

    if (std::holds_alternative<std::shared_ptr<DBNode>>(value)) {
        collect_entry_nodes_from_node(std::get<std::shared_ptr<DBNode>>(value), out);
        return;
    }
}

std::vector<std::shared_ptr<DBNode>> extract_nuisance_entries(const DBNode::Value& value) {
    std::vector<std::shared_ptr<DBNode>> out;
    collect_entry_nodes_from_value(value, out);
    return out;
}

}

NuisanceReader::NuisanceReader(std::shared_ptr<INuisancePathsProvider> paths_provider)
    : paths_provider_(std::move(paths_provider))
{
    if (!paths_provider_) {
        throw std::invalid_argument("NuisanceReader: paths_provider is null");
    }
}

fs::path NuisanceReader::default_path() const {
    return paths_provider_->default_nuisances_path();
}

fs::path NuisanceReader::user_path() const {
    return paths_provider_->user_nuisances_path();
}

NuisanceRegistry NuisanceReader::load_default() const {
    NuisanceRegistry registry;
    const fs::path path = default_path();

    if (path.empty()) {
        throw std::runtime_error("NuisanceReader: default path is empty");
    }
    if (!fs::exists(path)) {
        throw std::runtime_error("NuisanceReader: default file not found: " + path.string());
    }

    merge_file_into_registry(path, registry);
    return registry;
}

NuisanceRegistry NuisanceReader::load_user() const {
    return load_user(user_path());
}

NuisanceRegistry NuisanceReader::load_user(const fs::path& path) const {
    NuisanceRegistry registry;

    if (path.empty()) {
        return registry;
    }

    if (!fs::exists(path)) {
        throw std::runtime_error("NuisanceReader: user file not found: " + path.string());
    }

    merge_file_into_registry(path, registry);
    return registry;
}

NuisanceRegistry NuisanceReader::load() const {
    NuisanceRegistry registry = load_default();
    NuisanceRegistry user_registry = load_user();

    for (const auto& [pid, spec] : user_registry) {
        registry[pid] = spec;
    }

    return registry;
}


void NuisanceReader::merge_file_into_registry(const fs::path& path,
                                              NuisanceRegistry& registry) const
{
    auto provider = DBNodeProviderFactory::createDBNodeProvider(path);
    if (!provider) {
        throw std::runtime_error("NuisanceReader: could not create provider for: " + path.string());
    }

    auto root = provider->provide_db_as_node();
    if (!root) {
        throw std::runtime_error("NuisanceReader: provider returned null DBNode for: " + path.string());
    }

    merge_node_into_registry(*root, registry);
}

void NuisanceReader::merge_node_into_registry(const DBNode& root,
                                              NuisanceRegistry& registry) const
{
    if (!root.contains("nuisances")) {
        return;
    }

    const auto nuisances_value = root.get("nuisances");
    const auto entries = extract_nuisance_entries(nuisances_value);

    // if (entries.empty()) {
    //     std::ostringstream oss;
    //     oss << "NuisanceReader: 'nuisances' found but no entries could be extracted "
    //         << "(stored as " << value_kind(nuisances_value) << ")";
    //     throw std::runtime_error(oss.str());
    // }

    for (const auto& entry_ptr : entries) {
        if (!entry_ptr) {
            continue;
        }

        NuisanceSpec spec = parse_entry(*entry_ptr);
        registry[spec.param_id] = std::move(spec);
    }
}

NuisanceSpec NuisanceReader::parse_entry(const DBNode& entry) {
    NuisanceSpec spec;

    spec.param_id = make_param_id(entry);

    const double min_val = value_to_double(
        get_required_value(entry, {"min_val", "min"}),
        "min_val"
    );

    const double max_val = value_to_double(
        get_required_value(entry, {"max_val", "max"}),
        "max_val"
    );

    if (min_val > max_val) {
        std::ostringstream oss;
        oss << "NuisanceReader: invalid bounds for parameter (min_val="
            << min_val << " > max_val=" << max_val << ")";
        throw std::runtime_error(oss.str());
    }

    spec.bounds = {min_val, max_val};

    const std::string distribution = value_to_string(
        get_required_value(entry, {"distribution", "marginal"}),
        "distribution"
    );

    spec.marginal = parse_marginal_type(distribution);

    return spec;
}

ParamId NuisanceReader::make_param_id(const DBNode& entry) {
    const std::string block_str = value_to_string(
        get_required_value(entry, {"block"}),
        "block"
    );

    const std::string code_str = value_to_code_string(
        get_required_value(entry, {"code"}),
        "code"
    );

    BlockName block(block_str);
    LhaID code(code_str);

    return ParamId(block, code);
}

DBNode::Value NuisanceReader::get_required_value(const DBNode& node,
                                                 std::initializer_list<const char*> candidate_keys)
{
    for (const char* key : candidate_keys) {
        if (node.contains(key)) {
            return node.get(key);
        }
    }

    std::ostringstream oss;
    oss << "NuisanceReader: missing required field among {";
    bool first = true;
    for (const char* key : candidate_keys) {
        if (!first) oss << ", ";
        oss << key;
        first = false;
    }
    oss << "}";

    throw std::runtime_error(oss.str());
}

std::string NuisanceReader::value_to_string(const DBNode::Value& value,
                                            const std::string& field_name)
{
    if (std::holds_alternative<BlockName>(value)) {
        return strip_optional_quotes(std::get<BlockName>(value));
    }

    if (std::holds_alternative<int>(value)) {
        return std::to_string(std::get<int>(value));
    }

    if (std::holds_alternative<double>(value)) {
        std::ostringstream oss;
        oss << std::get<double>(value);
        return oss.str();
    }

    throw std::runtime_error(
        "NuisanceReader: field '" + field_name + "' must be string-like"
    );
}

std::string NuisanceReader::value_to_code_string(const DBNode::Value& value,
                                                 const std::string& field_name)
{
    const std::string s = value_to_string(value, field_name);
    return strip_optional_quotes(s);
}

double NuisanceReader::value_to_double(const DBNode::Value& value,
                                       const std::string& field_name)
{
    if (std::holds_alternative<double>(value)) {
        return std::get<double>(value);
    }

    if (std::holds_alternative<int>(value)) {
        return static_cast<double>(std::get<int>(value));
    }

    if (std::holds_alternative<BlockName>(value)) {
        const std::string s = strip_optional_quotes(std::get<BlockName>(value));
        if (is_numeric_string(s)) {
            return std::stod(s);
        }
    }

    throw std::runtime_error(
        "NuisanceReader: field '" + field_name + "' must be numeric"
    );
}

MarginalType NuisanceReader::parse_marginal_type(std::string raw) {
    raw = normalise(strip_optional_quotes(std::move(raw)));

    if (raw == "gaussian") {
        return MarginalType::GAUSSIAN;
    }

    if (raw == "half_gaussian" || raw == "halfgaussian" || raw == "split_gaussian") {
        return MarginalType::HALF_GAUSSIAN;
    }

    if (raw == "flat" || raw == "uniform") {
        return MarginalType::FLAT;
    }

    if (raw == "likelihood") {
        return MarginalType::LIKELIHOOD;
    }

    throw std::runtime_error("NuisanceReader: unknown distribution '" + raw + "'");
}

std::string NuisanceReader::normalise(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    std::replace(s.begin(), s.end(), '-', '_');
    std::replace(s.begin(), s.end(), ' ', '_');

    return s;
}