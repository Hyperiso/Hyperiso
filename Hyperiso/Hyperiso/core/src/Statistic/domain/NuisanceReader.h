#ifndef NUISANCEREADER_H
#define NUISANCEREADER_H

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_set>
#include <initializer_list>

#include "DBNode.h"
#include "DBNodeProviderFactory.h"
#include "INuisancePathsProvider.h"
#include "INuisanceReader.h"
#include "NuisanceSpec.h"

namespace fs = std::filesystem;

class NuisanceReader : public INuisanceReader {
public:
    explicit NuisanceReader(std::shared_ptr<INuisancePathsProvider> paths_provider);

    NuisanceRegistry load_default() const override;
    NuisanceRegistry load_user() const override;
    NuisanceRegistry load_user(const fs::path& user_path) const override;

    NuisanceRegistry load() const;

    fs::path default_path() const;
    fs::path user_path() const;

private:
    std::shared_ptr<INuisancePathsProvider> paths_provider_;

private:
    void merge_file_into_registry(const fs::path& path, NuisanceRegistry& registry) const;
    void merge_node_into_registry(const DBNode& root, NuisanceRegistry& registry) const;

    static NuisanceSpec parse_entry(const DBNode& entry);
    static ParamId make_param_id(const DBNode& entry);

    static DBNode::Value get_required_value(const DBNode& node,
                                            std::initializer_list<const char*> candidate_keys);

    static std::string value_to_string(const DBNode::Value& value,
                                       const std::string& field_name);

    static std::string value_to_code_string(const DBNode::Value& value,
                                            const std::string& field_name);

    static double value_to_double(const DBNode::Value& value,
                                  const std::string& field_name);

    static MarginalType parse_marginal_type(std::string raw);
    static std::string normalise(std::string s);
};

#endif // NUISANCEREADER_H