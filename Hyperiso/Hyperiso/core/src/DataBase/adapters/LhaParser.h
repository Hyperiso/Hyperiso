#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <string>
#include <map>
#include <memory>
#include <filesystem>
#include <climits>
#include <regex>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "IParser.h"
#include "lha_blocks.h"
#include "lha_elements.h"
#include "Logger.h"

/**
 * @enum TokenType
 * @brief Enumeration of token types used in parsing LHA files.
 */
enum class TokenType {
    FLOAT,  /**< Floating point number. */
    INTEGER,/**< Integer number. */
    BLOCK,  /**< Block keyword. */
    DECAY,  /**< Decay keyword. */
    NEWLINE,/**< Newline character. */
    SKIP,   /**< Whitespace character. */
    COMMENT,/**< Comment line. */
    WORD,   /**< Word token. */
    OTHER   /**< Other types of tokens. */
};

/**
 * @var analyzer_rx
 * @brief Regular expression for tokenizing LHA file content.
 */
const std::regex analyzer_rx(
    R"x(((?:[+-])?(?:\d+\.\d*|\.\d+)(?:[eEdD][+-]\d+)?|(\d+(?:[eEdD][+-]\d+)?))|(?:[+-]?\d+(?!\.))|(block)|(decay)|(\n)|([ \t]+)|(#.*)|([\w\=\.]+)|([^#]*))x",
    std::regex_constants::icase
); 

/**
 * @struct Token
 * @brief Represents a parsed token from the LHA file.
 */
struct Token {
    TokenType type;     /**< Type of the token. */
    std::string value;  /**< Text value of the token. */
    int row;            /**< Row number where the token is located. */
    int col;            /**< Column number where the token is located. */
};

/**
 * @class LhaParser
 * @brief Reads and manages LHA blocks and elements from LHA or FLHA files.
 */
class LhaParser : public IParser {
private:
    std::unordered_set<Prototype> blockPrototypes;                     /**< List of block prototypes used for parsing. */

    void fill_prototypes(std::string_view lha_path);

    /**
     * @brief Tokenizes the source string into individual tokens.
     */
    std::vector<Token> tokenize(const std::string& src) const;

    /**
     * @brief Parses the source string into blocks, optionally including comments.
     * @param comments If `true`, includes comments in the parsing process.
     */
    std::map<BlockName, std::vector<std::vector<std::string>>> parse_tokens(std::vector<Token> tokens, bool comments = false) const;

    /**
     * @brief Adds a new block to the reader from parsed lines.
     * @param id Block identifier.
     * @param lines Vector of lines containing the block's data.
     */
    void addBlock(std::map<BlockName, std::shared_ptr<LhaBlock>>& blocks, const std::string& id, const std::vector<std::vector<std::string>>& lines) const;

    /**
     * @brief Finds a prototype by block name.
     * @param name Name of the block.
     * @return `Prototype` instance if found; otherwise, an empty prototype.
     */
    Prototype findPrototype(std::string name) const;

    /**
     * @brief Converts the reader's content to a string representation.
     * @return String representing all blocks and their contents.
     */
    std::shared_ptr<Node> toDBNode(std::map<BlockName, std::shared_ptr<LhaBlock>> blocks) const;

public:
    /**
     * @brief Reads all blocks from the LHA file.
     */
    std::shared_ptr<Node> readFromFile(const std::string& input_file) const override;

    /**
     * @brief Reads all blocks from a string source.
     */
    std::shared_ptr<Node> parse(const std::string& src) const override;

    /**
     * @brief Writes the Node structure to a file.
     * @param filename The file name to write to.
     * @param root The root Node to write.
     */
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const;

    void set_prototypes(const std::unordered_set<Prototype>& prototypes);
};

#endif // HYPERISO_LHA_READER_H