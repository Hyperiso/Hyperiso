#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <fstream>
#include <string>
#include <map>
#include <memory>
#include <filesystem>

#include "lha_blocks.h"

class LhaReader {
private:
    std::map<BlockId, std::unique_ptr<LhaBlock>> blocks;
    bool isFLHA = false;
    std::filesystem::path lhaFile;

    void addBlock(BlockId id, std::ifstream& file);

public:
    LhaReader(std::string_view path);
    void readBlock(BlockId id);
    void readAll();

    LhaBlock* getBlock(BlockId id) {
        return blocks.at(id).get();
    }

    int getBlockCount() {
        return blocks.size();
    }
};


#endif // HYPERISO_LHA_READER_H