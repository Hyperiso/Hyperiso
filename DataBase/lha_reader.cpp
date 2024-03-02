#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

// #include "lha_blocks.h"
#include "lha_reader.h"

void Parser::tokenize() {
    int cLine = 0;
    int cCol = 0;

    auto rit = std::sregex_iterator(src.begin(), src.end(), analyzer_rx);
    auto rend = std::sregex_iterator();

    while (rit != rend) {
        std::smatch m = *rit;
        int group_index = m.size();
    
        for (int idx = 1; idx < m.size(); ++idx) {
            if (m[idx].matched) {
                if (idx == 1)
                    ++idx;
                group_index = idx - 2;
                break;
            }
        }

        auto tokenType = static_cast<TokenType>(group_index);
        auto value = m[group_index + 2].str();

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
    bool isWilson = false;
    bool isWilsonQ = false;
    std::string wilsonQ;
    std::string cBlock;

    int cCol = INT_MAX;

    for (const Token& t : this->tokens) {
        if (newBlock) {
            this->rawBlocks[t.value] = std::vector<std::vector<std::string>> {};
            cBlock = t.value;
            newBlock = false;
            isWilson = t.value.ends_with("WCOEF");
        } else if (t.type == TokenType::BLOCK) {
            newBlock = true;
        } else if (!comments && t.type == TokenType::COMMENT) {
            continue;
        } else if (isWilson && t.type == TokenType::WORD) {
            isWilsonQ = (t.value == "Q=");
        } else {
            if (isWilsonQ) {
                wilsonQ = t.value;
                isWilsonQ = false;
            }
            else {
                if (t.col < cCol) {
                    this->rawBlocks[cBlock].emplace_back(std::vector<std::string> {});
                    if(isWilson)
                        this->rawBlocks[cBlock].back().emplace_back(wilsonQ);
                }
                this->rawBlocks[cBlock].back().emplace_back(t.value);
                cCol = t.col;
            }
        }
    }
}

void LhaReader::addBlock(const std::string& id, const std::vector<std::vector<std::string>>& lines) {
    auto block = std::make_unique<LhaBlock>(id);
    block->readData(lines);
    this->blocks.insert(std::pair(id, std::move(block)));
}

LhaReader::LhaReader(std::string_view path) : lhaFile(std::filesystem::path(path)) {
    isFLHA = this->lhaFile.extension().string() == ".flha"; 
}

void LhaReader::readAll() {
    std::ifstream file(this->lhaFile.string());
    std::stringstream buffer;
    buffer << file.rdbuf();
    Parser parser {buffer.str()};
    parser.parse();
    auto blocks = parser.getBlocks();

    for (auto p : blocks) {
        addBlock(p.first, p.second);
    }
}

bool LhaReader::hasBlock(const std::string& id) const {
    return this->blocks.contains(id);  // C++20. Use blocks.count(id) != 0 before.
}

template <class T>
void LhaReader::extractFromBlock(std::string blockName, std::vector<T*>& vars) {
    LhaBlock* block = this->getBlock(blockName);
    if (block) {
        for (int id=0; id!=vars.size(); ++id) {
            *(vars.at(id)) = static_cast<LhaElement<T>*>(block->get(std::to_string(id + 1)))->getValue();
        }
    }
}

// int main() {

//     // std::ifstream file("../DataBase/example.flha");
//     // std::stringstream buffer;
//     // buffer << file.rdbuf();
//     // Parser parser {buffer.str()};
//     // parser.parse();
//     // auto blocks = parser.getBlocks();

//     // for (auto p : blocks) {
//     //     std::cout << "Block " << p.first << std::endl;
//     //     for (auto line: p.second) {
//     //         for (auto w: line) {
//     //             std::cout << w << '\t';
//     //         }
//     //         std::cout << '\n';
//     //     }
//     // }

//     LhaReader reader("../DataBase/example.flha");

//     reader.readAll();

//     std::cout << "Parsing ended, read " << reader.getBlockCount() << " block(s)." << std::endl;
//     // for (const auto& k : blockNames) {
//     //     BlockId id = BlockIdHelper::getBlockId(k);
//     //     if (reader.hasBlock(id)) {
//     //         std::cout << reader.getBlock(id)->toString() << std::endl;
//     //     }
//     // }

//     return 0;
// }
