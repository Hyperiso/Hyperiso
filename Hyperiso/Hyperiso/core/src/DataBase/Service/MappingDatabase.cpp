#include "MappingDatabase.h"
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"
#include <iostream>


MappingDatabase::MappingDatabase(const std::string& instanceName, const std::string& jsonFilePath,
                                 std::shared_ptr<IParamMappingSource> loader) {
    
    if (jsonFilePath.empty()) {
        std::cerr << "Instance " << instanceName
                  << " introuvable et aucun chemin JSON fourni." << std::endl;
    }
    
    if (!loader) {
        loader = std::make_shared<JsonParamMappingAdapter>(std::make_shared<JSONParser>());
    }

    load(jsonFilePath, loader);
}

void MappingDatabase::load(const std::string& jsonFilePath,
                           const std::shared_ptr<IParamMappingSource>& loader) {
    try {
        paramsMap = loader->loadFromFile(jsonFilePath);
    } catch (const std::exception& e) {
        std::cerr << "Echec chargement mapping depuis " << jsonFilePath
                  << " : " << e.what() << std::endl;
        paramsMap.clear();
    }
}

std::unordered_map<std::string, InterpretedParam> MappingDatabase::getParams() const {
    return std::unordered_map<std::string, InterpretedParam>(paramsMap.begin(), paramsMap.end());
}
