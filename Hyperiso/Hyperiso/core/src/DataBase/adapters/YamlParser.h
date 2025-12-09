#ifndef YAMLPARSER_H
#define YAMLPARSER_H

#include <memory>
#include <string>
#include <cctype>
#include <stdexcept>
#include <algorithm>
#include <stack>
#include "DBNode.h"
#include "IParser.h"
#include "Logger.h"

/**
 * @file YamlParser.h
 * @brief Lightweight YAML parser/serializer for Node trees.
 *
 * YAMLParser implements the IParser interface using a simple,
 * indentation-based subset of YAML. It supports:
 *   - nested mappings (key: value),
 *   - sequences introduced by '-' items,
 *   - scalar values (int, double, bool, string).
 *
 * The resulting structure is stored in the generic Node container.
 */

/**
 * @class YAMLParser
 * @brief YAML implementation of the IParser interface.
 *
 * YAMLParser can:
 *   - parse a YAML-like string into a Node hierarchy,
 *   - write a Node hierarchy to a YAML file,
 *   - read a Node hierarchy back from a YAML file.
 *
 * The parser is intentionally minimal and focuses on the subset of YAML
 * needed by the project (indentation-based maps and lists).
 */
class YAMLParser : public IParser {
public:
    /**
     * @brief Parses a YAML string into a Node hierarchy.
     *
     * @param input YAML text buffer.
     * @return Shared pointer to the root Node of the parsed document.
     *
     * @throws std::runtime_error on syntax or conversion errors.
     */
    std::shared_ptr<Node> parse(const std::string& input) const override;

    /**
     * @brief Writes a Node hierarchy to a YAML file.
     *
     * This uses Node::printYAMLToStream() to serialize the tree in a
     * readable YAML-like format.
     *
     * @param filename Target file name.
     * @param root     Root node to serialize.
     *
     * @throws std::runtime_error if the file cannot be opened or writing fails.
     */
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override;

    /**
     * @brief Reads and parses YAML data from a file.
     *
     * The file contents are loaded into memory and passed to parse().
     *
     * @param filename Source file name.
     * @return Shared pointer to the root Node of the parsed document.
     *
     * @throws std::runtime_error if the file cannot be opened or parsing fails.
     */
    std::shared_ptr<Node> readFromFile(const std::string& filename) const override;

private:
    /**
     * @brief Parses a single "key: value" YAML line.
     *
     * Splits @p line into a key and value part if a colon is present.
     * Leading and trailing spaces are stripped.
     *
     * @param line  Input line.
     * @param key   Output key (if any).
     * @param value Output value (if any).
     * @return true if a key/value could be extracted, false otherwise.
     */
    bool parseLine(const std::string& line, std::string& key, std::string& value) const;

    /**
     * @brief Adjusts node and indentation stacks when indentation changes.
     *
     * This helper is used by streaming parsers that track a stack of
     * active nodes and the associated indentation levels.
     *
     * @param indent      Current indentation of the line being parsed.
     * @param nodeStack   Stack of nodes corresponding to nested mappings.
     * @param indentStack Parallel stack of indentation levels.
     */
    void adjustIndentation(size_t indent, std::vector<std::shared_ptr<Node>>& nodeStack, std::vector<int>& indentStack) const;

    /**
     * @brief Recursively parses a YAML node from a stream.
     *
     * This function:
     *   - reads lines from @p stream,
     *   - uses indentation to determine when to return,
     *   - dispatches between key/value pairs and list items.
     *
     * @param stream      Input stream (already loaded with YAML content).
     * @param indentLevel Current indentation level for this node.
     * @return Shared pointer to the parsed Node.
     */
    static std::shared_ptr<Node> parseYAMLNode(std::istringstream& stream, int indentLevel);

    /**
     * @brief Processes a "key: value" line and attaches it to a Node.
     *
     * If the value is empty, a nested node is parsed at the next
     * indentation level. Otherwise, the scalar is parsed and stored
     * directly under the key.
     *
     * @param line        Current YAML line ("key: value" form).
     * @param node        Node to which the key/value is attached.
     * @param stream      Input stream for reading nested content.
     * @param indentLevel Indentation to use for nested nodes.
     */
    static void processKeyValue(const std::string& line, std::shared_ptr<Node>& node, std::istringstream& stream, int indentLevel);

    /**
     * @brief Processes a YAML list (sequence) under the current node.
     *
     * List items are introduced by '-'. Items can be:
     *   - scalars,
     *   - nested mappings (if followed by "key: value" structures).
     *
     * The resulting list is encoded as a group with numeric string keys
     * ("0", "1", ...) inside the target node.
     *
     * @param node        Node to which the list is attached.
     * @param stream      Input stream for reading further list items.
     * @param indentLevel Expected indentation for items.
     * @param firstLine   Optional first line already read (starting with '-').
     */
    static void processList(std::shared_ptr<Node>& node, std::istringstream& stream, int indentLevel, std::string firstLine = "");

    /**
     * @brief Trims leading and trailing spaces and tabs from a string.
     *
     * @param str Input string.
     * @return Trimmed string.
     */

    static std::string trim(const std::string& str);

    /**
     * @brief Counts leading spaces in a string.
     *
     * This is used to determine indentation level for YAML.
     *
     * @param str Input string.
     * @return Number of leading spaces (or npos if the line is all spaces).
     */
    static int countLeadingSpaces(const std::string& str);

    /**
     * @brief Detects whether a Node represents a "list-like" structure.
     *
     * A node is considered list-like if all of its keys are numeric
     * strings (e.g. "0", "1", "2"...). This is useful when deciding
     * how to render or interpret a Node.
     *
     * @param node Node to inspect.
     * @return true if the node looks like a list, false otherwise.
     */
    bool isListNode(const std::shared_ptr<Node>& node) const;

    /**
     * @brief Parses a numeric string as a double.
     *
     * @param value String to convert.
     * @return Parsed double.
     *
     * @throws std::runtime_error if the string cannot be converted.
     */
    double parseNumber(const std::string& value) const;

    /**
     * @brief Parses a scalar YAML value into a Node::Value.
     *
     * The order of detection is:
     *   - integer,
     *   - double,
     *   - boolean ("true"/"false"),
     *   - otherwise: string.
     *
     * @param value Raw scalar text.
     * @return Parsed Node::Value.
     */
    static Node::Value parseScalar(const std::string& value);

    /**
     * @brief Higher-level scalar parsing used in some code paths.
     *
     * Similar to parseScalar(), but uses slightly different detection
     * logic (e.g. looking for '.', 'e', 'E' to detect floating point).
     *
     * @param value Raw scalar text.
     * @return Parsed Node::Value.
     */
    Node::Value parseValue(const std::string& value) const;

    /**
     * @brief Checks if a string represents an integer.
     *
     * A string is considered integer if it is non-empty and all
     * characters are decimal digits.
     *
     * @param str Input string.
     * @return true if @p str is an integer, false otherwise.
     */
    static bool isDouble(const std::string& str);

    /**
     * @brief Checks if a string represents a floating-point number.
     *
     * Uses std::strtod under the hood to verify the entire string
     * can be parsed as a double.
     *
     * @param str Input string.
     * @return true if @p str is a valid double, false otherwise.
     */
    static bool isInteger(const std::string& str);
};

#endif // YAMLPARSER_H
