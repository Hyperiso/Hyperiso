#include "BlockAccessor.h"

void BlockAccessor::addBlock(const std::string& name, std::shared_ptr<Block> block) {
    blocks[name] = std::move(block);
}

double BlockAccessor::getValue(const std::string& blockName, int pdgCode) const {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getValue(pdgCode);
    }
    throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(pdgCode));
}

bool BlockAccessor::exist(const std::string blockName, int pdgCode) const {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return true;
    }
    throw false;
}

void BlockAccessor::setValue(const std::string& blockName, int pdgCode, double value, bool force) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        it->second->setValue(pdgCode, value, force);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

void BlockAccessor::setMode(const std::string& blockName, int pdgCode, ParameterMode mode) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        it->second->setMode(pdgCode, mode);
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::map<int, double> BlockAccessor::getAllValues(std::string blockName) {
    auto it = blocks.find(blockName);
    if (it != blocks.end()) {
        return it->second->getAllValues();
    } else {
        throw std::invalid_argument("Block not found");
    }
}

std::vector<std::string> BlockAccessor::get_blocks() {
    std::vector<std::string> keys;
    for (auto key : blocks){
        keys.push_back(key.first);
    }
    return keys;
}


