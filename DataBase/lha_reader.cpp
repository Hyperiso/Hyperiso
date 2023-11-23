#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <memory>
#include <vector>
#include <filesystem>
#include <algorithm>

#include "lha_blocks.h"
#include "lha_reader.h"

void LhaReader::addBlock(BlockId id, std::ifstream& file) {
    auto block = LhaBlockFactory::createBlock(id, isFLHA);
    if (block) {
        block->readData(file);
        blocks[id] = std::move(block);
    }
}

LhaReader::LhaReader(std::string_view path) : lhaFile(std::filesystem::path(path)) {
        isFLHA = this->lhaFile.extension().string() == ".flha"; 
    }

void LhaReader::readBlock(BlockId id) {
    if (this->blocks.contains(id)) return; // C++20. Use blocks.count(id) != 0 before.

    std::ifstream file(this->lhaFile.string());
    std::string line;
    std::string targetName = BlockIdHelper::getBlockName(id);

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

void LhaReader::readAll() {
    std::ifstream file(this->lhaFile.string());
    std::string line;
    while (std::getline(file, line)) {
        if (tolower(line[0]) != 'b') continue;
        std::string _, blockName;
        std::istringstream stream(line);
        stream >> _ >> blockName;
        addBlock(BlockIdHelper::getBlockId(blockName), file);
    }
}

int main() {
    LhaReader reader("../DataBase/example.flha");
    reader.readAll();

    std::cout << "Parsing ended, read " << reader.getBlockCount() << " block(s)." << std::endl;
    std::cout << reader.getBlock(BlockId::FCINFO)->toString() << std::endl;
    std::cout << reader.getBlock(BlockId::FMASS)->toString() << std::endl;

    return 0;
}
