#include "Json.h"

int main() {
    
    auto parser = ParserFactory::createParser(ParserFactory::Type::JSON);

    auto root = std::make_shared<Node>();

    root->set("wowow", "1", "2", "3");
    root->set("wowo2w", "1", "2", "4");
    root->set(42.5, "level1", "level2", "sdf");
    root->set("value", "level1", "level2", "ERER");

    std::cout << "Structure JSON générée :\n";
    root->printJSON();
    std::cout << "\n";

    try {
        auto value = std::get<BlockName>(root->get("1", "2", "3"));
        std::cout << "Valeur de [1, 2, 3] : " << value << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << "\n";
    }

    try {
        auto group = root->getGroup({"1", "2"});
        std::cout << "Groupe sous [1, 2] :\n";
        for (const auto& [key, val] : group) {
            std::cout << "- " << key << " : ";
            if (std::holds_alternative<BlockName>(val)) {
                std::cout << std::get<BlockName>(val);
            } else if (std::holds_alternative<int>(val)) {
                std::cout << std::get<int>(val);
            } else if (std::holds_alternative<double>(val)) {
                std::cout << std::get<double>(val);
            }
            std::cout << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << "\n";
    }

    std::map<std::string, Node::Value> newGroup = {
        {"newKey1", "newValue1"},
        {"newKey2", 12345},
        {"newKey3", false},
    };
    root->setGroup({"1", "2", "5"}, newGroup);

    std::cout << "Structure JSON après setGroup :\n";
    root->printJSON();
    std::cout << "\n";

    std::cout << "Structure YAML générée :\n";
    root->printYAML();
    parser->writeToFile("output.json", root);

    try {
        auto newRoot = parser->readFromFile("output.json");
        std::cout << "\nFichier JSON lu avec succès.\n";
        newRoot->printJSON();
    } catch (const std::exception& e) {
        std::cerr << "Erreur lors de la lecture : " << e.what() << "\n";
    }

    return 0;
}