#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <string>
#include <map>
#include <memory>
#include <filesystem>
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
    "((?:[\\+\\-])?((?:\\d+\\.\\d*|\\.\\d+)(?:[eEdD][\\+\\-]\\d+)?|\\d+(?:[eEdD][\\+\\-]\\d+))(?!\\.))|((?:[\\+\\-])?\\d+(?!\\.))|(block)|(\\n)|([ \\t]+)|(#.*)|([\\w\\=\\.]+)|([^#]*)",
    std::regex_constants::icase
); 

struct Token {
    TokenType type;
    std::string value;
    int row;
    int col;
};

class Parser {
    public:
        inline explicit Parser(std::string src) : src(std::move(src)) {}
        void parse(bool comments = false);
        inline std::vector<std::vector<std::string>> getBlock(const std::string& blockName) { return this->rawBlocks[blockName]; }
        inline std::map<std::string, std::vector<std::vector<std::string>>> getBlocks() { return this->rawBlocks; }

    private:
        const std::string src;
        std::vector<Token> tokens;
        std::map<std::string, std::vector<std::vector<std::string>>> rawBlocks;

        void tokenize();
};

class LhaReader {
private:
    std::map<std::string, std::unique_ptr<LhaBlock>> blocks;
    bool isFLHA = false;
    std::filesystem::path lhaFile;

    void addBlock(const std::string& id, const std::vector<std::vector<std::string>>& lines);

public:
    LhaReader(std::string_view path);
    bool hasBlock(const std::string& id) const;
    void readAll();

    inline LhaBlock* getBlock(const std::string& id) const {
        return blocks.at(id).get();
    }

    inline int getBlockCount() const {
        return blocks.size();
    }
};


#endif // HYPERISO_LHA_READER_H