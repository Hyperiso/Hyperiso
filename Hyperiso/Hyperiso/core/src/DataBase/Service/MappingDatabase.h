#ifndef MAPPING_DATABASE_H
#define MAPPING_DATABASE_H

#include <unordered_map>
#include <string>
#include <memory>
#include "IParamMappingSource.h" 

class MappingDatabase {
public:

    MappingDatabase(const std::string& instanceName, const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader);

    std::unordered_map<std::string, InterpretedParam> getParams() const;

private:

    void load(const std::string& jsonFilePath,
              const std::shared_ptr<IParamMappingSource>& loader);

private:
    std::unordered_map<std::string, InterpretedParam> paramsMap;
};

#endif // MAPPING_DATABASE_H
