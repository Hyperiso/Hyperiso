#ifndef JSONPARAMMAPPINGADAPTER_H
#define JSONPARAMMAPPINGADAPTER_H

#include "IParamMappingSource.h"
#include "DBNode.h"
#include <variant>
#include <iostream>

class JsonParamMappingAdapter : public IParamMappingSource {
public:
    explicit JsonParamMappingAdapter(std::shared_ptr<IParser> parser)
        : parser_(std::move(parser)) {}

    std::unordered_map<std::string, InterpretedParam>
    loadFromFile(const std::string& filePath) const override;

private:
    static std::string asString(const Node::Value& v) {
        if (auto s = std::get_if<BlockName>(&v)) return *s;
        if (auto i = std::get_if<int>(&v))      return std::to_string(*i);
        if (auto d = std::get_if<double>(&v))   return std::to_string(*d);
        if (auto b = std::get_if<bool>(&v))     return *b ? "true" : "false";
        throw std::runtime_error("JSON mapping: valeur non scalaire convertissable en string");
    }

private:
    std::shared_ptr<IParser> parser_;
};

#endif // JSONPARAMMAPPINGADAPTER_H
