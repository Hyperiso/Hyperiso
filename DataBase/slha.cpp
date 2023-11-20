#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <memory>

// Interface de base pour les blocs
class SlhaBlock {
public:
    virtual void readData(std::ifstream& file) = 0;
    // Autres méthodes communes
};

// Exemple de bloc spécifique (le plus simple)
class MassBlock : public SlhaBlock {
    // Données spécifiques au bloc MASS
public:
    void readData(std::ifstream& file) override {
        // Lire les données du bloc MASS
    }
};

// Factory pour créer des instances de blocs
class SlhaBlockFactory {
public:
    static std::unique_ptr<SlhaBlock> createBlock(const std::string& blockName, bool isFLHA) {
        std::string effectiveBlockName = isFLHA ? "F" + blockName : blockName;
        if (effectiveBlockName == "MASS" || effectiveBlockName == "FMASS") {
            return std::make_unique<MassBlock>();
        }
        //  On ajoute les blocs ici
        return nullptr;
    }
};

// Classe principale pour lire les fichiers SLHA
class ReadSlha {
    std::map<std::string, std::vector<std::unique_ptr<SlhaBlock>>> blocks;
    bool isFLHA = false;

public:
    void readFromFile(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        // Détection du type de fichier (SLHA ou FLHA)
        // ...

        while (std::getline(file, line)) {
            std::string blockName = /* Analyse de la ligne */;
            auto block = SlhaBlockFactory::createBlock(blockName, isFLHA);
            if (block) {
                block->readData(file);
                blocks[blockName].push_back(std::move(block));
            }
        }
    }

    // Méthodes pour accéder aux données des blocs
};

int main() {
    ReadSlha reader;
    reader.readFromFile("example.slha");

    return 0;
}
