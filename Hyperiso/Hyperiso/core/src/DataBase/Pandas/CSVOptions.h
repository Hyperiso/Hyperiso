
#ifndef CSVOPTIONS_H
#define CSVOPTIONS_H
#include <string>
#include <vector>
#include <typeindex>
#include <unordered_map>

struct CSVOptions {
    bool hasIndex = false;
    std::unordered_map<std::string, std::type_index> columnTypes;
    std::type_index indexType = typeid(std::string);

    std::type_index& operator[](const std::string& colName) {
        auto result = columnTypes.emplace(colName, std::type_index(typeid(void)));
        return result.first->second;
    }
};

#endif