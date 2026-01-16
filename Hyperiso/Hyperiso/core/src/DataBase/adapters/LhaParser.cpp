#include "LhaParser.h"

void LhaParser::addBlock(std::map<BlockName, std::shared_ptr<LhaBlock>>& blocks, const BlockName& id, const std::vector<std::vector<std::string>>& lines) const {
    auto block = std::make_shared<LhaBlock>(findPrototype(id));
    LOG_DEBUG(id);
    block->readData(lines);
    BlockName id_ci = block->getPrototype().blockName;
    id_ci.to_upper();
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


std::map<BlockName, std::vector<std::vector<std::string>>>
LhaParser::parse_tokens(std::vector<Token> tokens, bool comments) const
{
    bool newBlock       = false;
    bool hasGlobalScale = false;
    bool isQ            = false;
    bool skipBlock      = false;
    bool decay          = false;

    std::string globalQ;
    BlockName   cBlock;
    int         cCol = INT_MAX;

    std::map<BlockName, std::vector<std::vector<std::string>>> rawBlocks;
    auto current_block_it = rawBlocks.end();

    for (const Token& t : tokens) {
        if (newBlock) {
            auto prototype = this->findPrototype(t.value);
            if (prototype.blockName != "") {
                LOG_DEBUG("LHA reader: Block ", prototype.blockName, " found.");

                cBlock        = prototype.blockName;
                hasGlobalScale= prototype.globalScale;
                skipBlock     = false;
                decay         = false;
                globalQ.clear();
                cCol = INT_MAX;

                auto [it, inserted] =
                    rawBlocks.emplace(cBlock, std::vector<std::vector<std::string>>{});
                current_block_it = it;
            }
            else if (decay) {
                LOG_DEBUG("LHA reader: Decay block found. Skipping.");
                skipBlock        = true;
                hasGlobalScale   = false;
                current_block_it = rawBlocks.end();
                decay            = false;
            }
            else {
                LOG_WARN("LHA reader: Unknown block " + t.value + " encountered. Skipping.");
                skipBlock        = true;
                hasGlobalScale   = false;
                current_block_it = rawBlocks.end();
            }
            newBlock = false;
        }
        else if (t.type == TokenType::BLOCK) {
            newBlock = true;
        }
        else if (t.type == TokenType::DECAY) {
            decay    = true;
            newBlock = true;
        }
        else if (!comments && t.type == TokenType::COMMENT) {
            continue;
        }
        else if (hasGlobalScale && t.type == TokenType::WORD && !skipBlock) {
            isQ = (t.value == "Q=" || t.value == "q=");
        }
        else if (!skipBlock && current_block_it != rawBlocks.end()) {
            if (isQ) {
                globalQ = t.value;
                isQ     = false;
            }
            else {
                LOG_VERBOSE("Token : [", static_cast<int>(t.type), ", ", t.value, "]");
                auto& rows = current_block_it->second;

                if (t.col <= cCol) {
                    rows.emplace_back(std::vector<std::string>{});
                    if (hasGlobalScale)
                        rows.back().emplace_back(!globalQ.empty() ? globalQ : "-1");
                }
                rows.back().emplace_back(t.value);
                cCol = t.col;
            }
        }
    }

    return rawBlocks;
}

Prototype LhaParser::findPrototype(BlockName name) const {
    for (const auto& p : blockPrototypes) {
        name.to_upper();
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
    std::map<BlockName, std::shared_ptr<LhaBlock>> blocks;
    for (auto &[id, lines] : rawBlocks) {
        addBlock(blocks, id, lines);
    }     
    LOG_DEBUG("LHA file parsed.");

    return this->toDBNode(blocks);
}

void LhaParser::writeToFile(const std::string &filename,
                            const std::shared_ptr<Node> &root) const {
    // TODO   
    std::string a = filename;
    a+std::string("a");
    root->get("truc");                             
}

void LhaParser::set_prototypes(const std::unordered_set<Prototype> &prototypes)
{
    this->blockPrototypes = prototypes;
}


std::shared_ptr<Node> LhaParser::toDBNode(std::map<BlockName, std::shared_ptr<LhaBlock>> blocks) const {
    Node root;

    for (const auto& [blockName, blockPtr] : blocks) {
        auto block_node = blockPtr->toDBNode();

        auto group = block_node->getGroup({blockName});
        root.setGroup({blockName}, group);

        if (block_node->contains("scale")) {
            auto sval = block_node->get("scale");
            if (std::holds_alternative<double>(sval)) {
                root.set(std::get<double>(sval), blockName, "scale");
            } else if (std::holds_alternative<int>(sval)) {
                root.set(static_cast<double>(std::get<int>(sval)), blockName, "scale");
            } else {
                LOG_WARN("Expected numeric 'scale' for block ", blockName, " but got variant index ", sval.index());
            }
        }
    }

    return std::make_shared<Node>(root);
}
