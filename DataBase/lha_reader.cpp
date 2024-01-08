#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <memory>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <regex>

#include "lha_blocks.h"
#include "lha_reader.h"

std::vector<Token> Tokenizer::tokenize() {
    auto rit = std::sregex_iterator(src.begin(), src.end(), analyzer_rx);
    auto rend = std::sregex_iterator();

    while (rit != rend) {
        
    } 
}

std::unique_ptr<LhaBlock> LhaReader::addBlock(BlockId id, std::ifstream& file) {
    auto block = LhaBlockFactory::createBlock(id, isFLHA);
    if (block) {
        block->readData(file);
    }
    return block;
}

LhaReader::LhaReader(std::string_view path) : lhaFile(std::filesystem::path(path)) {
    isFLHA = this->lhaFile.extension().string() == ".flha"; 
}

std::unique_ptr<LhaBlock> LhaReader::readBlock(BlockId id) {
    if (this->hasBlock(id)) return nullptr; 

    std::ifstream file(this->lhaFile.string());
    std::string line;
    std::string targetName = BlockIdHelper::getBlockName(id);

    while (std::getline(file, line)) {
        if (tolower(line[0]) != 'b') continue;
        std::string _, name;
        std::istringstream stream(line);
        stream >> _ >> name;
        if (name == targetName) {
            return addBlock(id, file);
        }
    }
    return nullptr;
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

bool LhaReader::hasBlock(BlockId id) const {
    return this->blocks.contains(id);  // C++20. Use blocks.count(id) != 0 before.
}

int main() {
    LhaReader reader("../DataBase/example.flha");
    reader.readAll();

    std::cout << "Parsing ended, read " << reader.getBlockCount() << " block(s)." << std::endl;
    for (const auto& k : blockNames) {
        BlockId id = BlockIdHelper::getBlockId(k);
        if (reader.hasBlock(id)) {
            std::cout << reader.getBlock(id)->toString() << std::endl;
        }
    }

    return 0;
}
