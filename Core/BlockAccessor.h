#pragma once

#include <map>
#include <string>
#include <memory>
#include <stdexcept>
#include "Block.h"

enum class FlavorParamType {
    LIFETIME,
    DECAY_CONSTANT,
    DECAY_CONSTANT_RATIO
};

class BlockAccessor : public Block {
public:
    void addBlock(const std::string& name, std::shared_ptr<Block> block) {
        blocks[name] = std::move(block);
    }

    double getValue(const std::string& blockName, int pdgCode) const {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            return it->second->getValue(pdgCode);
        }
        throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(pdgCode));
    }
    
    bool exist(const std::string blockName, int pdgCode) const {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            return true;
        }
        throw false;
    }

    void setValue(const std::string& blockName, int pdgCode, double value, bool force = false) {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            it->second->setValue(pdgCode, value, force);
        } else {
            throw std::invalid_argument("Block not found");
        }
    }

    void setMode(const std::string& blockName, int pdgCode, ParameterMode mode) {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            it->second->setMode(pdgCode, mode);
        } else {
            throw std::invalid_argument("Block not found");
        }
    }

    std::map<int, double> getAllValues(std::string blockName) {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            return it->second->getAllValues();
        } else {
            throw std::invalid_argument("Block not found");
        }
    }

    std::vector<std::string> get_blocks() {
        std::vector<std::string> keys;
        for (auto key : blocks){
            keys.push_back(key.first);
        }
        return keys;
    }
    std::map<int, double> getAllValues() override {
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    // These methods are to satisfy the Block interface, but can be unused directly
    double getValue(int pdgCode) const override {
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    void setValue(int pdgCode, double value, bool force = false) override {
        throw std::logic_error("Use setValue with block name for BlockAccessor");
    }

    void setMode(int pdgCode, ParameterMode mode) override {
        throw std::logic_error("Use setMode with block name for BlockAccessor");
    }

private:
    std::map<std::string, std::shared_ptr<Block>> blocks;
};