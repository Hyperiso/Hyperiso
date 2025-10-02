// MappingDatabase.cpp
#include "MappingDatabase.h"
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"      // ton JSONParser concret
#include <iostream>

std::unordered_map<std::string, std::shared_ptr<MappingDatabase>> MappingDatabase::instances;

std::shared_ptr<MappingDatabase> MappingDatabase::getInstance(
    const std::string& instanceName,
    const std::string& jsonFilePath,
    std::shared_ptr<IParamMappingSource> loader
) {
    auto it = instances.find(instanceName);
    if (it != instances.end()) return it->second;

    if (jsonFilePath.empty()) {
        std::cerr << "Instance " << instanceName
                  << " introuvable et aucun chemin JSON fourni." << std::endl;
        return nullptr;
    }

    // Adaptateur par défaut : JSON
    if (!loader) {
        loader = std::make_shared<JsonParamMappingAdapter>(std::make_shared<JSONParser>());
    }

    auto db = std::shared_ptr<MappingDatabase>(new MappingDatabase(jsonFilePath, loader));
    instances[instanceName] = db;
    return db;
}

MappingDatabase::MappingDatabase(const std::string& jsonFilePath,
                                 std::shared_ptr<IParamMappingSource> loader) {
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

const std::unordered_map<std::string, InterpretedParam>& MappingDatabase::getParams() const {
    return paramsMap;
}
