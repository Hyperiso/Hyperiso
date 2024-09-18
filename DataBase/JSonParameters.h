#ifndef JSON_PARAM_PARSER_H
#define JSON_PARAM_PARSER_H
#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

class JSONBlock {
public:
    void setValue(int pdgCode, double value) {
        pdgCodeValues[pdgCode] = value;
    }

    double getValue(int pdgCode) const {
        auto it = pdgCodeValues.find(pdgCode);
        if (it != pdgCodeValues.end()) {
            return it->second;
        } else {
            throw std::runtime_error("pdgCode non trouvé dans le block.");
        }
    }

    std::string toJSON(int indentLevel = 1) const;

    void fromJSON(const std::string& json);

private:
    std::unordered_map<int, double> pdgCodeValues;
};

class JSONParser {
public:
    void addElement(const std::string& blockName, int pdgCode, double value) {
        blocks[blockName].setValue(pdgCode, value);
    }

    double getElement(const std::string& blockName, int pdgCode) const {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            return it->second.getValue(pdgCode);
        } else {
            throw std::runtime_error("Block non trouvé.");
        }
    }

    std::string toJSON(int indentLevel = 0) const;

    void fromJSON(const std::string& json);

    void saveToFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Impossible d'ouvrir le fichier pour écrire.");
        }
        file << toJSON();
        file.close();
    }

    void loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Impossible d'ouvrir le fichier pour lire.");
        }
        std::stringstream buffer;
        buffer << file.rdbuf();
        fromJSON(buffer.str());
        file.close();
    }

    static JSONParser* getInstance(int id) {
        if (instances.find(id) == instances.end()) {
            instances[id] = new JSONParser();
        }
        return instances[id];
    }

    static void removeInstance(int id) {
    // Supprime l'instance avec l'identifiant donné
    instances.erase(id);
    }

private:
    JSONParser() = default;

    static std::unordered_map<int, JSONParser*> instances;
    std::unordered_map<std::string, JSONBlock> blocks;
};

#endif