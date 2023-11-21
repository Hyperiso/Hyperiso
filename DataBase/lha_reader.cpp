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
    std::vector<std::unique_ptr<AbstractElement>> entries;
    BlockId id;

public:
    LhaBlock(BlockId id) : id(id) {}

    virtual void readData(std::ifstream& file) = 0;

    AbstractElement* get(std::string_view id) const { 
        auto p = [id](const std::unique_ptr<AbstractElement>& e) { 
            return e->getId() == id; 
        };
        return (*std::find_if(entries.begin(), entries.end(), p)).get(); 
    }

    std::string toString() const {
        std::stringstream stream;
        stream << "Block " << getBlockName(this->id) << ":\n";
        for (auto& entry: entries) {
            stream << entry->toString();
        }
        return stream.str();
    }
};

class MassBlock : public LhaBlock {
public:
    MassBlock(BlockId id) : LhaBlock(id) {}

    void readData(std::ifstream& file) override {
        std::string line, id;
        int scheme;
        double value, scale;
        while (std::getline(file, line)) {
            if (line.at(0) == '#') continue;
            std::istringstream stream(line);
            stream >> id >> value >> scheme >> scale;
            RenormalizationScheme s = static_cast<RenormalizationScheme>(scheme);
            this->entries.emplace_back(std::make_unique<GeneralElement<double>>(id, value, s, scale));
            if (tolower(file.peek()) == 'b' || file.peek() == EOF) break;
        }
    }
};

class InfoBlock : public LhaBlock {
public:
    InfoBlock(BlockId id) : LhaBlock(id) {}

    void readData(std::ifstream& file) override {
        std::string line, id, value;
        while (std::getline(file, line)) {
            if (line.at(0) == '#') continue;
            std::istringstream stream(line);
            stream >> id >> value;
            this->entries.emplace_back(std::make_unique<TypedElement<std::string>>(id, value));
            if (tolower(file.peek()) == 'b' || file.peek() == EOF) break;
        }
    }
}; 

// class ModSelBlock : public LhaBlock {
// public:
//     ModSelBlock(BlockId id) : LhaBlock(id) {}

//     void readData(std::ifstream& file) override {
//         std::string line;
//         while (std::getline(file, line)) {
//             if (line.at(0) == '#') continue;
//             std::istringstream stream(line);
//             int id;
//             if (id == 99) {
//                 std::string value;
//                 stream >> value;
//                 this->entries.emplace_back(std::make_unique<TypedElement<std::string>>(id, value));
//             }
            
//             if (tolower(file.peek()) == 'b' || file.peek() == EOF) break;
//         }
//     }
// }; 

// Factory pour créer des instances de blocs
class LhaBlockFactory {
public:
    static std::unique_ptr<LhaBlock> createBlock(BlockId id, bool isFLHA) {
        if (id == BlockId::MASS && !isFLHA || id == BlockId::FMASS && isFLHA) {
            return std::make_unique<MassBlock>(id);
        } else if (id == BlockId::FCINFO && isFLHA) {
            return std::make_unique<InfoBlock>(id);
        }
        // TODO : Parse all block names (switch ?)
        return nullptr;
    }
};


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
        isFLHA = this->lhaFile.extension().string() == ".flha"; 
    }

    void readBlock(BlockId id) {
        if (this->blocks.contains(id)) return; // C++20. Use blocks.count(id) != 0 before.

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

    LhaBlock* getBlock(BlockId id) {
        return blocks.at(id).get();
    }

    int getBlockCount() {
        return blocks.size();
    }
};

int main() {
    LhaReader reader("../DataBase/example.flha");
    reader.readAll();

    std::cout << "Parsing ended, read " << reader.getBlockCount() << " block(s)." << std::endl;
    std::cout << reader.getBlock(BlockId::FCINFO)->toString() << std::endl;
    std::cout << reader.getBlock(BlockId::FMASS)->toString() << std::endl;

    return 0;
}
