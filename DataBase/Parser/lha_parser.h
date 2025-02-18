#ifndef __LHA_PARSER_H__
#define __LHA_PARSER_H__

#include <string>
#include <map>
#include <memory>
#include <filesystem>
#include <regex>
#include <climits>

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
 * @class Parser
 * @brief Parses LHA files and extracts blocks of data.
 */
class Parser {
    public:
        /**
         * @brief Constructs a Parser with the given source and reader.
         * @param src Source string to be parsed.
         */
        inline explicit Parser(std::string src) : src(std::move(src)) {}

        /**
         * @brief Retrieves parsed tokens.
         * @return Vector of tokens read in the source file.
         */
        inline std::vector<Token> getTokens() { return this->tokens; }

        /**
         * @brief Tokenizes the source string into individual tokens.
         */
        void tokenize();


    private:
        const std::string src;                                                  /**< Source string to be parsed. */
        std::vector<Token> tokens;                                              /**< Tokens parsed from the source string. */
};

#endif // __LHA_PARSER_H__
