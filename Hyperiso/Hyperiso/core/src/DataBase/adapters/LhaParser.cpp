#include "LhaParser.h"

void LhaParser::addBlock(std::map<std::string, std::shared_ptr<LhaBlock>>& blocks, const std::string& id, const std::vector<std::vector<std::string>>& lines) const {
    auto block = std::make_shared<LhaBlock>(findPrototype(id));
    LOG_DEBUG(id);
    block->readData(lines);
    std::string id_ci = id;
    std::transform(id_ci.begin(), id_ci.end(), id_ci.begin(), ::toupper);
    blocks.insert(std::pair(id_ci, std::move(block)));
}

std::vector<Token> LhaParser::tokenize(const std::string &src) const {
    int cLine = 0;
    int cCol = 0;

    auto rit = std::sregex_iterator(src.begin(), src.end(), analyzer_rx);
    auto rend = std::sregex_iterator();

    std::vector<Token> tokens;
    while (rit != rend) {
        std::smatch m = *rit;
        size_t group_index = m.size();

        for (size_t idx = 1; idx < m.size(); ++idx) {
            if (m[idx].matched) {
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
            tokens.emplace_back(Token{tokenType, value, cLine, cCol});
            ++cCol;
        }

        ++rit;
    }

    return tokens;
}

std::map<std::string, std::vector<std::vector<std::string>>> LhaParser::parse_tokens(std::vector<Token> tokens, bool comments) const {
    bool newBlock = false;
    bool hasGlobalScale = false;
    bool isQ = false;
    bool skipBlock = false;
    bool decay = false;
    std::string globalQ;
    std::string cBlock;

    int cCol = INT_MAX;

    std::map<std::string, std::vector<std::vector<std::string>>> rawBlocks;

    for (const Token& t : tokens) {
        if (newBlock) {
            auto prototype = this->findPrototype(t.value);
            if (prototype.blockName != "") {
                LOG_DEBUG("LHA reader: Block " + prototype.blockName + " found.");
                rawBlocks[t.value] = std::vector<std::vector<std::string>> {};
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
                    rawBlocks[cBlock].emplace_back(std::vector<std::string> {});
                    if(hasGlobalScale)
                        rawBlocks[cBlock].back().emplace_back(globalQ != "" ? globalQ : "-1");
                }
                rawBlocks[cBlock].back().emplace_back(t.value);
                cCol = t.col;
            }
        }
    }

    return rawBlocks;
}

Prototype LhaParser::findPrototype(std::string name) const {
    for (const auto& p : blockPrototypes) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (p.blockName == name) 
            return p;
    }
    return Prototype{""};
}

std::shared_ptr<Node> LhaParser::readFromFile(const std::string& input_file) const {
    std::ifstream file(input_file.data());
    std::stringstream buffer;
    buffer << file.rdbuf();
    return this->parse(buffer.str());
}

std::shared_ptr<Node> LhaParser::parse(const std::string &src) const {
    auto rawBlocks = this->parse_tokens(this->tokenize(src));
    std::map<std::string, std::shared_ptr<LhaBlock>> blocks;
    for (auto &[id, lines] : rawBlocks) {
        addBlock(blocks, id, lines);
    }     
    LOG_DEBUG("LHA file parsed.");
    return this->toDBNode(blocks);
}

void LhaParser::writeToFile(const std::string &filename,
                            const std::shared_ptr<Node> &root) const {
                                
}

void LhaParser::set_prototypes(const std::unordered_set<Prototype> &prototypes)
{
    this->blockPrototypes = prototypes;
}

std::shared_ptr<Node> LhaParser::toDBNode(std::map<std::string, std::shared_ptr<LhaBlock>> blocks) const {
    Node node;
    for (const auto& block : blocks) {
        node.set(block.second->toDBNode(), "lha_root");
    }
    return std::make_shared<Node>(node);
}
