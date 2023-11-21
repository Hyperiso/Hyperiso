#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <memory>
#include <vector>
#include <filesystem>
#include <algorithm>

#include "lha_elements.cpp"

// TODO : Add SLHA block names 
enum class BlockId {  FCINFO, FMODSEL, SMINPUTS, VCKMIN, UPMNSIN, VCKM, IMVCKM,
                        UPMNS, IMUPMNS, FMASS, MASS, FLIFE, FCONST, FCONSTRATIO, FBAG,
                        FWCOEF, IMFWCOEF, FOBS, FOBSERR, FOBSSM, FDIPOLE, FPARAM};

const std::vector<std::string> blockNames = {
    "FCINFO", "FMODSEL", "SMINPUTS", "VCKMIN", "UPMNSIN",
    "VCKM", "IMVCKM", "UPMNS", "IMUPMNS", "FMASS", "MASS",
    "FLIFE", "FCONST", "FCONSTRATIO", "FBAG", "FWCOEF",
    "IMFWCOEF", "FOBS", "FOBSERR", "FOBSSM", "FDIPOLE", "FPARAM"
};

std::string getBlockName(BlockId value) {
    return blockNames[static_cast<int>(value)];
}

BlockId getBlockId(const std::string& value) {
    int idx = std::find(blockNames.begin(), blockNames.end(), value) - blockNames.begin();
    return static_cast<BlockId>(idx);
}


// Interface de base pour les blocs
class LhaBlock {
protected:
    std::map<int, LhaElement> content;   

public:
    virtual void readData(std::ifstream& file) = 0;
    LhaElement get(int id) { return content.at(id); }
};

// Exemple de bloc spécifique (le plus simple)
class MassBlock : public LhaBlock {
public:
    void readData(std::ifstream& file) override {
        std::string line;
        while (std::getline(file, line)) {
            if (tolower(file.peek()) == 'b' || file.peek() == EOF) break;
            int id;
            double value;
            std::istringstream stream(line);
            stream >> id >> value;
            LhaReal elt (id, value);
            this->content.insert({id, elt});
        }
    }
};

// Factory pour créer des instances de blocs
class LhaBlockFactory {
public:
    static std::unique_ptr<LhaBlock> createBlock(BlockId id, bool isFLHA) {
        if (id == BlockId::MASS && !isFLHA || id == BlockId::FMASS && isFLHA) {
            return std::make_unique<MassBlock>();
        }
        //  On ajoute les blocs ici
        return nullptr;
    }
};

// Classe principale pour lire les fichiers SLHA
class LhaReader {
private:
    std::map<BlockId, std::unique_ptr<LhaBlock>> blocks;
    bool isFLHA = false;
    std::filesystem::path lhaFile;

    void addBlock(BlockId id, std::ifstream& file) {
        auto block = LhaBlockFactory::createBlock(id, isFLHA);
        if (block) {
            block->readData(file);
            blocks[id] = std::move(block);
        }
    }

public:
    LhaReader(std::string_view path) : lhaFile(std::filesystem::path(path)) {
        if (this->lhaFile.extension().string() == ".flha") {
            isFLHA = true;
        }
    }

    void readBlock(BlockId id) {
        std::ifstream file(this->lhaFile.string());
        std::string line;
        std::string targetName = getBlockName(id);

        while (std::getline(file, line)) {
            if (tolower(line[0]) != 'b') continue;
            std::string _, name;
            std::istringstream stream(line);
            stream >> _ >> name;
            if (name == targetName) {
                addBlock(id, file);
                break;
            }
        }
    }

    void readAll() {
        std::ifstream file(this->lhaFile.string());
        std::string line;
        while (std::getline(file, line)) {
            if (tolower(line[0]) != 'b') continue;
            std::string _, blockName;
            std::istringstream stream(line);
            stream >> _ >> blockName;
            addBlock(getBlockId(blockName), file);
        }
    }

    // Méthodes pour accéder aux données des blocs

    LhaBlock* getBlock(BlockId id) {
        return blocks.at(id).get();
    }

    int getBlockCount() {
        return blocks.size();
    }
};

int main() {
    LhaReader reader("../DataBase/example.flha");
    reader.readBlock(BlockId::FMASS);

    std::cout << "Parsing ended, read " << reader.getBlockCount() << " block(s)." << std::endl;
    std::cout << reader.getBlock(BlockId::FMASS)->get(13).getValue() << std::endl;

    return 0;
}
