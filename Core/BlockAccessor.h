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

// Composite pattern to manage access to multiple blocks
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
        std::cout<< blockName << std::endl;
        throw std::invalid_argument("Block " + blockName + " not found with pdg code : " + std::to_string(pdgCode));
    }

    void setValue(const std::string& blockName, int pdgCode, double value) {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            it->second->setValue(pdgCode, value);
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

    // These methods are to satisfy the Block interface, but can be unused directly
    double getValue(int pdgCode) const override {
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    void setValue(int pdgCode, double value) override {
        throw std::logic_error("Use setValue with block name for BlockAccessor");
    }

    void setMode(int pdgCode, ParameterMode mode) override {
        throw std::logic_error("Use setMode with block name for BlockAccessor");
    }

private:
    std::map<std::string, std::shared_ptr<Block>> blocks;
};

class FlavorBlockAccessor : public FlavorBlock {
public:
    void addBlock(FlavorParamType name, std::shared_ptr<FlavorBlock> block) {
        blocks[name] = std::move(block);
    }

    double getValue(FlavorParamType blockName, std::string id) const {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            return it->second->getValue(id);
        }
        throw std::invalid_argument("Block not found");
    }

    void setValue(FlavorParamType blockName, std::string id, double value) {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            it->second->setValue(id, value);
        } else {
            throw std::invalid_argument("Block not found");
        }
    }

    void setMode(FlavorParamType blockName, std::string id, ParameterMode mode) {
        auto it = blocks.find(blockName);
        if (it != blocks.end()) {
            it->second->setMode(id, mode);
        } else {
            throw std::invalid_argument("Block not found");
        }
    }
private:
    std::map<FlavorParamType, std::shared_ptr<FlavorBlock>> blocks;
};