#ifndef YAML_INPUT_READER_H
#define YAML_INPUT_READER_H

#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <variant>

#include "YamlParser.h"
#include "DBNode.h"
#include "LhaID.h"

struct YamlInputParam {
    std::string block_name;
    LhaID pdg_code;
};

struct YamlScanParam {
    std::string block_name;
    LhaID pdg_code;
    double min_val;
    double max_val;
    double step_val;
};

class YamlInputReader {
public:
    explicit YamlInputReader(const std::string& filename);

    std::vector<YamlInputParam> get_input_params() const;
    std::vector<YamlScanParam>  get_scan_params() const;

private:
    YAMLParser yaml_p;
    std::shared_ptr<DBNode> root_;

    static std::string unquote(std::string s);

    static std::string as_string(const DBNode::Value& v, const std::string& ctx);
    static double      as_double(const DBNode::Value& v, const std::string& ctx);
    static std::shared_ptr<DBNode> as_node(const DBNode::Value& v, const std::string& ctx);
};

#endif