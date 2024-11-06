#ifndef MAPPING_DATABASE_H
#define MAPPING_DATABASE_H

#include <unordered_map>
#include <string>
#include <memory>

struct InterpretedParam {
    std::string block;
    int pdgCode;
};

class MappingDatabase {
public:
    static std::shared_ptr<MappingDatabase> getInstance(const std::string& instanceName, const std::string& jsonFilePath = "");

    const std::unordered_map<std::string, InterpretedParam>& getParams() const;

private:
    MappingDatabase(const std::string& jsonFilePath);
    void loadFromJson(const std::string& jsonFilePath);

    std::unordered_map<std::string, InterpretedParam> paramsMap;
    static std::unordered_map<std::string, std::shared_ptr<MappingDatabase>> instances;
};

#endif // MAPPING_DATABASE_H
