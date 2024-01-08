#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <fstream>
#include <string>
#include <map>
#include <memory>
#include <filesystem>
#include <optional>
#include <regex>

#include "lha_blocks.h"

enum class TokenType {
    FLOAT,
    INTEGER,
    BLOCK,
    NEWLINE,
    SKIP,
    COMMENT,
    WORD,
    OTHER
};

const std::regex analyzer_rx(
    "((?:[\\+\\-])?((?:\\d+\\.\\d*|\\.\\d+)(?:[eEdD][\\+\\-]\\d+)?|\\d+(?:[eEdD][\\+\\-]\\d+))(?!\\.))|((?:[\\+\\-])?\\d+(?!\\.))|(^[a-z]+)|(\\n)|([ \\t]+)|(#.)|([\\w\\=\\.]+)|([^#]*)",
    std::regex_constants::icase
); 

struct Token {
    TokenType type;
    std::optional<std::string> value;
};

class Tokenizer {
    public:
        inline explicit Tokenizer(std::string src) : src(std::move(src)) {}
        std::vector<Token> tokenize();

    private:
        const std::string src;
        int index;

        std::optional<TokenType> peek() const;
        Token consume();
};

class LhaReader {
private:
    std::map<BlockId, std::unique_ptr<LhaBlock>> blocks;
    bool isFLHA = false;
    std::filesystem::path lhaFile;

    std::unique_ptr<LhaBlock> addBlock(BlockId id, std::ifstream& file);

public:
    LhaReader(std::string_view path);
    bool hasBlock(BlockId id) const;
    std::unique_ptr<LhaBlock> readBlock(BlockId id);
    void readAll();

    inline LhaBlock* getBlock(BlockId id) const {
        return blocks.at(id).get();
    }

    inline int getBlockCount() const {
        return blocks.size();
    }
};


#endif // HYPERISO_LHA_READER_H