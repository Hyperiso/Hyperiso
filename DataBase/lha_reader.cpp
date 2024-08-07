#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "lha_reader.h"
#include "Logger.h"

void Parser::tokenize() {
    int cLine = 0;
    int cCol = 0;

    auto rit = std::sregex_iterator(src.begin(), src.end(), analyzer_rx);
    auto rend = std::sregex_iterator();

    while (rit != rend) {
        std::smatch m = *rit;
        size_t group_index = m.size();

        // Don't touch, it works.
        for (size_t idx = 1; idx < m.size(); ++idx) {
            if (m[idx].matched) {
                // if (idx == 1)  
                //     ++idx;
                group_index = idx - 1;
                break;
            }
        }

        auto tokenType = static_cast<TokenType>(group_index);
        auto value = m[group_index + 1].str();

        if (tokenType == TokenType::NEWLINE) {
            ++cLine;
            cCol = 0;
        } else if (tokenType != TokenType::SKIP && value != "") {
            this->tokens.emplace_back(Token{tokenType, value, cLine, cCol});
            ++cCol;
        }

        ++rit;
    }
}

void Parser::parse(bool comments) {
    this->tokenize();
    bool newBlock = false;
    bool hasGlobalScale = false;
    bool isQ = false;
    bool skipBlock = false;
    bool decay = false;
    std::string globalQ;
    std::string cBlock;

    int cCol = INT_MAX;

    for (const Token& t : this->tokens) {
        if (newBlock) {
            auto prototype = reader->findPrototype(t.value);
            if (prototype.blockName != "") {
                LOG_INFO("LHA reader: Block " + prototype.blockName + " found.");
                this->rawBlocks[t.value] = std::vector<std::vector<std::string>> {};
                cBlock = t.value;
                hasGlobalScale = prototype.globalScale;
                skipBlock = false;
            } else if (decay) {
                LOG_INFO("LHA reader: Decay block found. Skipping.");
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
                // std::cout << "Token : [" << (int)t.type << ", " << t.value  << "]" << std::endl;
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

void LhaReader::addBlock(const std::string& id, const std::vector<std::vector<std::string>>& lines) {
    auto block = std::make_unique<LhaBlock>(findPrototype(id));
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

void LhaReader::readAll() {
    std::ifstream file(this->lhaFile.string());
    std::stringstream buffer;
    buffer << file.rdbuf();
    Parser parser {buffer.str(), this};
    parser.parse();
    auto blocks = parser.getBlocks();
    
    for (auto p : blocks) {
        addBlock(p.first, p.second);
    }
          
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
    LOG_INFO("Updating LHA blocks...");
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
