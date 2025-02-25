#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "lha_reader.h"
#include "Logger.h"
#include "lha_parser.h"

void LhaReader::addBlock(const std::string& id, const std::vector<std::vector<std::string>>& lines) {
    auto block = std::make_shared<LhaBlock>(findPrototype(id));
    LOG_DEBUG(id);
    block->readData(lines);
    std::string id_ci = id;
    std::transform(id_ci.begin(), id_ci.end(), id_ci.begin(), ::toupper);
    this->blocks.insert(std::pair(id_ci, std::move(block)));
}

LhaReader::LhaReader(std::string_view path) : lhaFile(std::filesystem::path(path)) {
    isFLHA = this->lhaFile.extension().string() == ".flha"; 
    blockPrototypes = SLHA_BLOCKS;
    if (isFLHA) {
        blockPrototypes.insert(blockPrototypes.end(), FLHA_BLOCKS.begin(), FLHA_BLOCKS.end());
    }
}

void LhaReader::parse_tokens(std::vector<Token> tokens, bool comments) {
    bool newBlock = false;
    bool hasGlobalScale = false;
    bool isQ = false;
    bool skipBlock = false;
    bool decay = false;
    std::string globalQ;
    std::string cBlock;

    int cCol = INT_MAX;

    for (const Token& t : tokens) {
        if (newBlock) {
            auto prototype = this->findPrototype(t.value);
            if (prototype.blockName != "") {
                LOG_DEBUG("LHA reader: Block " + prototype.blockName + " found.");
                this->rawBlocks[t.value] = std::vector<std::vector<std::string>> {};
                cBlock = t.value;
                hasGlobalScale = prototype.globalScale;
                skipBlock = false;
            } else if (decay) {
                LOG_DEBUG("LHA reader: Decay block found. Skipping.");
                skipBlock = true;
                decay = false;
            } else {
                LOG_WARN("LHA reader: Unknown block " + t.value + " encountered. Skipping.");
                skipBlock = true;
            }
            newBlock = false;
        } else if (t.type == TokenType::BLOCK) {
            newBlock = true;
        } else if (t.type == TokenType::DECAY) {
            decay = true;
            newBlock = true;
        } else if (!comments && t.type == TokenType::COMMENT) {
            continue;
        } else if (hasGlobalScale && t.type == TokenType::WORD && !skipBlock) {
            isQ = (t.value == "Q=" || t.value == "q=");
        } else if (!skipBlock) {
            if (isQ) {
                globalQ = t.value;
                isQ = false;
            }
            else {
                LOG_VERBOSE("Token : [", (int)t.type, ", ", t.value, "]");
                if (t.col <= cCol) {   
                    this->rawBlocks[cBlock].emplace_back(std::vector<std::string> {});
                    if(hasGlobalScale)
                        this->rawBlocks[cBlock].back().emplace_back(globalQ != "" ? globalQ : "-1");
                }
                this->rawBlocks[cBlock].back().emplace_back(t.value);
                cCol = t.col;
            }
        }
    }
}

bool LhaReader::hasElement(const std::string &block_id, const LhaID &elt_id) const {
    return this->hasBlock(block_id) && this->getBlock(block_id)->hasElement(elt_id);
}

void LhaReader::readAll() {
    std::ifstream file(this->lhaFile.string());
    std::stringstream buffer;
    buffer << file.rdbuf();
    LhaParser parser {buffer.str()};
    parser.tokenize();
    this->parse_tokens(parser.getTokens());
    
    for (auto &[id, lines] : this->rawBlocks) {
        addBlock(id, lines);
    }
          
    LOG_DEBUG("LHA file parsed.");
}

Prototype LhaReader::findPrototype(std::string name) const {
    for (auto p : blockPrototypes) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (p.blockName == name) 
            return p;
    }
    return Prototype{""};
}

std::string LhaReader::getLhaPath() const {
    return this->lhaFile.string();
}

void LhaReader::update(std::string_view newLha) {
    LOG_DEBUG("Updating LHA blocks...");
    this->lhaFile = std::filesystem::path(newLha);
    isFLHA = lhaFile.extension().string() == ".flha";
    
    if (isFLHA && this->blockPrototypes.size() == SLHA_BLOCKS.size()) {
        this->blockPrototypes.insert(blockPrototypes.end(), FLHA_BLOCKS.begin(), FLHA_BLOCKS.end());
    }
    this->blocks.clear();
     
    this->readAll();
}

std::string LhaReader::toString() const
{
    std::stringstream ss;
    for (const auto& block: blocks) {
        ss << "\n--- Block ID: " << block.first << " ---\n";
        ss << block.second->toString();
    }

    return ss.str();
}

bool LhaReader::hasBlock(const std::string& id) const {
    return this->blocks.contains(id);  // C++20. Use blocks.count(id) != 0 before.
}
