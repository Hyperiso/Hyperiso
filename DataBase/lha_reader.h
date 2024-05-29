#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <string>
#include <map>
#include <memory>
#include <filesystem>
#include <regex>

#include "lha_blocks.h"
#include "lha_elements.h"

enum class TokenType {
    FLOAT,
    INTEGER,
    BLOCK,
    DECAY,
    NEWLINE,
    SKIP,
    COMMENT,
    WORD,
    OTHER
};

const std::regex analyzer_rx(
    R"x(((?:[+-])?(?:\d+\.\d*|\.\d+)(?:[eEdD][+-]\d+)?|(\d+(?:[eEdD][+-]\d+)?))|(?:[+-]?\d+(?!\.))|(block)|(decay)|(\n)|([ \t]+)|(#.*)|([\w\=\.]+)|([^#]*))x",
    std::regex_constants::icase
); 

struct Token {
    TokenType type;
    std::string value;
    int row;
    int col;
};

class LhaReader;

class Parser {
    public:
        inline explicit Parser(std::string src, LhaReader* reader) : src(std::move(src)), reader(reader) {}
        void parse(bool comments = false);
        inline std::vector<std::vector<std::string>> getBlock(const std::string& blockName) { return this->rawBlocks[blockName]; }
        inline std::map<std::string, std::vector<std::vector<std::string>>> getBlocks() { return this->rawBlocks; }

    private:
        const std::string src;
        const LhaReader* reader;
        std::vector<Token> tokens;
        std::map<std::string, std::vector<std::vector<std::string>>> rawBlocks;


        void tokenize();
};

class LhaReader {
private:
    std::vector<Prototype> blockPrototypes;
    std::map<std::string, std::unique_ptr<LhaBlock>> blocks;
    bool isFLHA = false;
    std::filesystem::path lhaFile;

    void addBlock(const std::string& id, const std::vector<std::vector<std::string>>& lines);

public:
    LhaReader(std::string_view path);
    bool hasBlock(const std::string& id) const;
    void readAll();
    Prototype findPrototype(std::string name) const;
    std::string getLhaPath() const;
    void update(std::string_view newLha);

    inline void addBlockType(std::string blockName, int itemCount=2, int valueIdx=1, int scaleIdx=-1, int rgIdx=-1, bool globalScale=false) {
        std::transform(blockName.begin(), blockName.end(), blockName.begin(), ::toupper);  // Make sure block name is uppercase 
        this->blockPrototypes.emplace_back(Prototype{blockName, itemCount, valueIdx, scaleIdx, rgIdx, globalScale});
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T*>& vars) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (int id=0; id!=vars.size(); ++id) {
                *(vars.at(id)) = static_cast<LhaElement<T>*>(block->get(std::to_string(id + 1)))->getValue();
            }
        }
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T>& vars) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (int id=0; id < vars.size(); ++id) {
                vars.at(id) = static_cast<LhaElement<T>*>(block->get(std::to_string(id + 1)))->getValue();
            }
        }
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T>& vars, const std::vector<int>& ids) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (int i=0; i < vars.size(); ++i) {
                auto e = block->get(std::to_string(ids.at(i)));
                vars.at(i) = e ? static_cast<LhaElement<T>*>(e)->getValue() : T {};
            }
        }
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T>& vars, std::vector<std::string>& ids) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (int i=0; i < vars.size(); ++i) {
                auto e = block->get(ids[i]);
                vars[i] = e ? static_cast<LhaElement<T>*>(e)->getValue() : T {};
            }
        }
    }

    inline LhaBlock* getBlock(std::string id) const {
        std::transform(id.begin(), id.end(), id.begin(), ::toupper);
        return this->hasBlock(id) ? blocks.at(id).get() : nullptr;
    }

    template <typename T>
    inline T getValue(const std::string& blockName, const std::string& eltId) {
        return static_cast<LhaElement<T>*>(this->getBlock(blockName)->get(eltId))->getValue();
    }

    inline int getBlockCount() const {
        return blocks.size();
    }

    std::string toString() const;
};


#endif // HYPERISO_LHA_READER_H