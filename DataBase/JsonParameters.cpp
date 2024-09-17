#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

class Block {
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

    std::string toJSON() const {
        std::stringstream ss;
        ss << "{";
        bool first = true;
        for (const auto& [pdgCode, value] : pdgCodeValues) {
            if (!first) ss << ",";
            ss << "\"" << pdgCode << "\":" << value;
            first = false;
        }
        ss << "}";
        return ss.str();
    }

    void fromJSON(const std::string& json) {
        pdgCodeValues.clear();
        std::istringstream stream(json);
        char c;
        stream >> c;
        while (stream >> c && c != '}') {
            if (c == '"') {
                int pdgCode;
                stream >> pdgCode; 
                stream >> c >> c;
                double value;
                stream >> value;
                pdgCodeValues[pdgCode] = value;
                stream >> c;
            }
        }
    }

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

    std::string toJSON() const {
        std::stringstream ss;
        ss << "{";
        bool first = true;
        for (const auto& [blockName, block] : blocks) {
            if (!first) ss << ",";
            ss << "\"" << blockName << "\":" << block.toJSON();
            first = false;
        }
        ss << "}";
        return ss.str();
    }

    void fromJSON(const std::string& json) {
        blocks.clear();
        std::istringstream stream(json);
        char c;
        stream >> c;
        while (stream >> c && c != '}') {
            if (c == '"') {
                std::string blockName;
                std::getline(stream, blockName, '"');
                stream >> c >> c;
                std::string blockData;
                std::getline(stream, blockData, '}');
                blockData += '}';
                Block block;
                block.fromJSON(blockData);
                blocks[blockName] = block;
                stream >> c;
            }
        }
    }

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

private:
    std::unordered_map<std::string, Block> blocks;
};

int main() {
    JSONParser parser;


    parser.addElement("block1", 2112, 10.5);
    parser.addElement("block1", 2212, 20.3);
    parser.addElement("block2", 11, 0.511);

    parser.saveToFile("data.json");

    JSONParser newParser;
    newParser.loadFromFile("data.json");

    try {
        std::cout << "block1, 2112 : " << newParser.getElement("block1", 2112) << std::endl;
        std::cout << "block1, 2212 : " << newParser.getElement("block1", 2212) << std::endl;
        std::cout << "block2, 11 : " << newParser.getElement("block2", 11) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }

    return 0;
}
