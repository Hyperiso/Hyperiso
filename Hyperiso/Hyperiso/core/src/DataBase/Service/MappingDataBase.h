#ifndef MAPPING_DATABASE_H
#define MAPPING_DATABASE_H

#include <unordered_map>
#include <string>
#include <memory>

#include "IParser.h"

/**
 * @struct InterpretedParam
 * @brief Represents a parameter with a block name and a PDG code.
 */
struct InterpretedParam {
    std::string block; /**< The block name associated with the parameter. */
    int pdgCode; /**< The PDG code associated with the parameter. */
};

/**
 * @class MappingDatabase
 */
class MappingDatabase {
public:
    /**
     * @brief Retrieves an instance of MappingDatabase by name, loading it from a JSON file if necessary.
     * 
     * Returns an existing instance of the database associated with the given `instanceName`,
     * or creates and loads a new instance from the specified JSON file if it doesn't already exist.
     *
     * @param instanceName The name of the instance.
     * @param jsonFilePath Optional file path to the JSON file used for loading the database.
     * @return Shared pointer to the MappingDatabase instance.
     */
    static std::shared_ptr<MappingDatabase> getInstance(const std::string& instanceName, const std::string& jsonFilePath = "");

    /**
     * @brief Retrieves the map of interpreted parameters.
     * 
     * Provides access to the internal mapping of parameter names to `InterpretedParam` objects.
     * 
     * @return Constant reference to the unordered map of parameters.
     */
    const std::unordered_map<std::string, InterpretedParam>& getParams() const;

private:
    /**
     * @brief Constructs the MappingDatabase and loads data from a JSON file.
     * 
     * Initializes the database and loads parameters from the specified JSON file.
     *
     * @param jsonFilePath Path to the JSON file used to initialize the database.
     */
    MappingDatabase(const std::string& jsonFilePath);

    /**
     * @brief Loads parameters from a JSON file into the database.
     * 
     * Parses the JSON file specified by `jsonFilePath` and populates the `paramsMap` with
     * extracted parameters and associated data.
     *
     * @param jsonFilePath Path to the JSON file.
     */
    void loadFromJson(const std::string& jsonFilePath);

    std::unordered_map<std::string, InterpretedParam> paramsMap; /**< Map of parameter names to `InterpretedParam` data. */
    static std::unordered_map<std::string, std::shared_ptr<MappingDatabase>> instances; /**< Map of instance names to MappingDatabase pointers. */
    
};

#endif // MAPPING_DATABASE_H
