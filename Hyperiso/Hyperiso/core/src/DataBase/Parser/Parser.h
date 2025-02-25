#ifndef PARSER_H
#define PARSER_H

#include <memory>
#include <string>
#include <cctype>
#include <stdexcept>
#include <algorithm>
#include <stack>
#include "DBNode.h"
#include "Logger.h"

/**
 * @brief Abstract base class for data parsers.
 */
class Parser {
public:
    virtual ~Parser() = default;

    /**
     * @brief Parses input into a Node structure.
     * @param input The string input to parse.
     * @return A shared pointer to the root Node.
     */
    virtual std::shared_ptr<Node> parse(const std::string& input) const = 0;

    /**
     * @brief Writes the Node structure to a file.
     * @param filename The file name to write to.
     * @param root The root Node to write.
     */
    virtual void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const = 0;

    /**
     * @brief Reads a Node structure from a file.
     * @param filename The file name to read from.
     * @return A shared pointer to the root Node.
     */
    virtual std::shared_ptr<Node> readFromFile(const std::string& filename) const = 0;
};

/**
 * @brief JSON implementation of the Parser class.
 */
class JSONParser : public Parser {
public:
    std::shared_ptr<Node> parse(const std::string& input) const override;
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override;
    std::shared_ptr<Node> readFromFile(const std::string& filename) const override;

    std::shared_ptr<Node> parseObject(const std::string& input, size_t& index) const;

    std::string parseString(const std::string& input, size_t& index) const;

    Node::Value parseValue(const std::string& input, size_t& index) const;

    double parseNumber(const std::string& input, size_t& index) const;

    void skipWhitespace(const std::string& input, size_t& index) const;

    std::shared_ptr<Node> parseArray(const std::string& input, size_t& index) const;

};

/**
 * @brief YAML implementation of the Parser class.
 */
class YAMLParser : public Parser {
public:
    std::shared_ptr<Node> parse(const std::string& input) const override;
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override;
    std::shared_ptr<Node> readFromFile(const std::string& filename) const override;

private:
    size_t countLeadingSpaces(const std::string& line) const;
    bool parseLine(const std::string& line, std::string& key, std::string& value) const;
    void adjustIndentation(size_t indent, std::vector<std::shared_ptr<Node>>& nodeStack, std::vector<int>& indentStack) const;
    void trim(std::string& str) const;
    Node::Value parseValue(const std::string& value) const;
};

/**
 * @brief Factory for creating Parser instances.
 */
class ParserFactory {
public:
    enum class Type { JSON, YAML };

    /**
     * @brief Creates a parser instance of the specified type.
     * @param type The type of parser to create.
     * @return A unique pointer to the created parser.
     */
    static std::unique_ptr<Parser> createParser(Type type);
};

#endif // PARSER_H