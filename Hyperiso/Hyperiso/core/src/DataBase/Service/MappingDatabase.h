// MappingDatabase.h
#ifndef MAPPING_DATABASE_H
#define MAPPING_DATABASE_H

#include <unordered_map>
#include <string>
#include <memory>
#include "IParamMappingSource.h" // définit InterpretedParam + port

class MappingDatabase {
public:
    // getInstance avec port/adaptateur injectable (par défaut: JSON)
    static std::shared_ptr<MappingDatabase> getInstance(
        const std::string& instanceName,
        const std::string& jsonFilePath = "",
        std::shared_ptr<IParamMappingSource> loader = nullptr
    );

    const std::unordered_map<std::string, InterpretedParam>& getParams() const;

private:
    MappingDatabase(const std::string& jsonFilePath,
                    std::shared_ptr<IParamMappingSource> loader);

    void load(const std::string& jsonFilePath,
              const std::shared_ptr<IParamMappingSource>& loader);

private:
    std::unordered_map<std::string, InterpretedParam> paramsMap;
    static std::unordered_map<std::string, std::shared_ptr<MappingDatabase>> instances;
};

#endif // MAPPING_DATABASE_H
