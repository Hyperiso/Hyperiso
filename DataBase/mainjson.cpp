#include "Json.h"

int main() {
    
    auto parser = ParserFactory::createParser(ParserFactory::Type::JSON);

    auto root = std::make_shared<Node>();
    root->set("wowow", "1", "2", "3");
    // root->set(42.5, "level1", "level2", "sdf");
    // root->set("value", "level1", "level2", "ERER");
    std::cout << "Structure JSON générée :\n";
    root->printJSON();
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