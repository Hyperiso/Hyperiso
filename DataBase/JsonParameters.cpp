#include "JSonParameters.h"

std::string JSONBlock::toJSON(int indentLevel) const {
    std::stringstream ss;
    std::string indent(indentLevel * 2, ' ');  // Indentation de 2 espaces par niveau
    ss << "{\n";
    bool first = true;
    for (const auto& [pdgCode, value] : pdgCodeValues) {
        if (!first) ss << ",\n";
        ss << indent << "\"" << pdgCode << "\": " << value;
        first = false;
    }
    ss << "\n" << std::string((indentLevel - 1) * 2, ' ') << "}"; // Fermeture du bloc avec indentation
    return ss.str();
}

void JSONBlock::fromJSON(const std::string& json) {
    pdgCodeValues.clear();
    std::istringstream stream(json);
    char c;
    int pdgCode;
    double value;
    while (stream >> c) {
        if (c == '"') {
            stream >> pdgCode;
            stream >> c >> c;
            stream >> value;
            pdgCodeValues[pdgCode] = value;
            stream >> c;
            if (c == '}') break;
        }
    }
}

std::string JSONParser::toJSON(int indentLevel) const {
    std::stringstream ss;
    std::string indent(indentLevel * 2, ' ');  // Indentation de 2 espaces par niveau
    ss << "{\n";
    bool first = true;
    for (const auto& [blockName, block] : blocks) {
        if (!first) ss << ",\n";
        ss << indent << "\"" << blockName << "\": " << block.toJSON(indentLevel + 1);
        first = false;
    }
    ss << "\n" << std::string(indentLevel * 2, ' ') << "}"; // Fermeture avec indentation
    return ss.str();
}

void JSONParser::fromJSON(const std::string& json) {
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
            JSONBlock block;
            block.fromJSON(blockData);
            blocks[blockName] = block;
            stream >> c;
        }
    }
}

std::unordered_map<int, JSONParser*> JSONParser::instances;

// int main() {
//     JSONParser parser;


//     parser.addElement("block1", 2112, 0.5);
//     parser.addElement("block1", 2212, 0.3);
//     parser.addElement("block2", 11, 0.5121);

//     parser.saveToFile("data.json");

//     JSONParser newParser;
//     newParser.loadFromFile("data.json");

//     try {
//         std::cout << "block1, 2112 : " << newParser.getElement("block1", 2112) << std::endl;
//         std::cout << "block1, 2212 : " << newParser.getElement("block1", 2212) << std::endl;
//         std::cout << "block2, 11 : " << newParser.getElement("block2", 11) << std::endl;
//         parser.addElement("block1", 211432, 0.5);
//         parser.saveToFile("data.json");
//     } catch (const std::exception& e) {
//         std::cerr << "Erreur : " << e.what() << std::endl;
//     }

//     return 0;
// }

// int main() {
//     // Création de deux instances de JSONParser avec des identifiants spécifiques
//     JSONParser* parser1 = JSONParser::getInstance(100);  // Instance avec ID 100
//     JSONParser* parser2 = JSONParser::getInstance(200);  // Instance avec ID 200

//     // Ajout d'éléments dans la première instance
//     parser1->addElement("block1", 2112, 10.5);
//     parser1->addElement("block1", 2212, 20.3);

//     // Ajout d'éléments dans la deuxième instance
//     parser2->addElement("block2", 11, 0.511);
//     parser2->addElement("block2", 13, 105.7);

//     // Sauvegarde des deux instances dans des fichiers distincts
//     parser1->saveToFile("data_parser100.json");
//     parser2->saveToFile("data_parser200.json");

//     // Chargement des données à partir des fichiers
//     JSONParser* loadedParser1 = JSONParser::getInstance(100);
//     loadedParser1->loadFromFile("data_parser100.json");

//     JSONParser* loadedParser2 = JSONParser::getInstance(200);
//     loadedParser2->loadFromFile("data_parser200.json");

//     // Affichage des éléments de la première instance chargée
//     std::cout << "Parser 100 - block1, 2112 : " << loadedParser1->getElement("block1", 2112) << std::endl;

//     // Affichage des éléments de la deuxième instance chargée
//     std::cout << "Parser 200 - block2, 11 : " << loadedParser2->getElement("block2", 11) << std::endl;

//     // Suppression de l'instance avec ID 100
//     JSONParser::removeInstance(100);

//     // Tentative d'accès à l'instance supprimée (ID 100)
//     try {
//         JSONParser* removedParser = JSONParser::getInstance(100);
//         std::cout << removedParser->getElement("block1", 2112) << std::endl;
//     } catch (const std::exception& e) {
//         std::cerr << "Erreur : " << e.what() << std::endl;
//     }

//     return 0;
// }
