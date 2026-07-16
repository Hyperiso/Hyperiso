#ifndef JSONPARSER_H
#define JSONPARSER_H

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
 * @file JsonParser.h
 * @brief Lightweight JSON parser/serializer for DBNode trees.
 *
 * JSONParser implements the IParser interface for a small, well-defined
 * subset of JSON:
 *   - objects: { "key": value, ... }
 *   - arrays:  [ value, value, ... ]
 *   - scalars: strings, numbers, booleans (true/false)
 *
 * Parsed data is stored in the generic DBNode hierarchy, using:
 *   - keys mapped directly to DBNode children,
 *   - arrays encoded as nodes with numeric string keys ("0", "1", ...).
 */

/**
 * @class JSONParser
 * @brief JSON implementation of the IParser interface.
 *
 * JSONParser can:
 *   - parse a JSON string into a DBNode tree,
 *   - serialize a DBNode tree to JSON and write it to a file (via DBNode::printJSONToStream),
 *   - read JSON from a file and convert it back into a DBNode tree.
 *
 * The parser is intentionally simple and does not support the full JSON
 * standard (e.g. no null literal, no escape sequences in strings).
 */
class JSONParser : public IParser {
public:
    /**
     * @brief Parses a JSON document from an in-memory string.
     *
     * The top-level JSON value is expected to be an object (`{ ... }`).
     *
     * @param input JSON text buffer.
     * @return Shared pointer to the root DBNode representing the JSON object.
     *
     * @throws std::runtime_error on malformed JSON or conversion errors.
     */
    std::shared_ptr<DBNode> parse(const std::string& input) const override;

    /**
     * @brief Serializes a DBNode tree to a JSON file.
     *
     * Internally uses DBNode::printJSONToStream() to produce a JSON-like
     * representation and writes it to @p filename.
     *
     * @param filename Target file name.
     * @param root     Root DBNode to serialize.
     *
     * @throws std::runtime_error if the file cannot be opened or a write error occurs.
     */
    void writeToFile(const std::string& filename, const std::shared_ptr<DBNode>& root) const override;

    /**
     * @brief Reads a JSON file and parses it into a DBNode tree.
     *
     * The entire file content is loaded into memory and passed to parse().
     *
     * @param filename Source file name.
     * @return Shared pointer to the root DBNode representing the JSON document.
     *
     * @throws std::runtime_error if the file cannot be opened or parsing fails.
     */
    std::shared_ptr<DBNode> readFromFile(const std::string& filename) const override;

    /**
     * @brief Parses a JSON object starting at the given index.
     *
     * Expects input[index] to point to a '{', and consumes the entire
     * object, returning a DBNode where each key becomes a child entry.
     *
     * Example JSON:
     *   { "a": 1, "b": true }
     *
     * becomes a DBNode with:
     *   - key "a" -> 1
     *   - key "b" -> true
     *
     * @param input JSON string.
     * @param index Reference to the current position in the string;
     *              advanced past the closing '}' on success.
     * @return Shared pointer to the DBNode representing the object.
     *
     * @throws std::runtime_error on invalid syntax.
     */
    std::shared_ptr<DBNode> parseObject(const std::string& input, size_t& index) const;

    /**
     * @brief Parses a JSON string literal.
     *
     * Expects input[index] to point to the opening double quote (`"`).
     * Reads until the next unescaped `"`. Escape sequences are not
     * currently interpreted; characters are copied verbatim.
     *
     * @param input JSON string.
     * @param index Reference to the current position; on return, points
     *              just after the closing `"`.
     * @return The extracted string value (without quotes).
     *
     * @throws std::runtime_error if the expected quotes are missing.
     */
    std::string parseString(const std::string& input, size_t& index) const;

    /**
     * @brief Parses any JSON value at the current position.
     *
     * Supported value types:
     *   - object:   { ... }
     *   - array:    [ ... ]
     *   - string:   "..."
     *   - number:   integer or floating-point
     *   - boolean:  true / false
     *
     * The result is converted into a DBNode::Value:
     *   - objects and arrays become DBNode subtrees,
     *   - scalars become int/double/bool/string stored in the variant.
     *
     * @param input JSON string.
     * @param index Reference to the current position; advanced past
     *              the parsed value on success.
     * @return Parsed DBNode::Value.
     *
     * @throws std::runtime_error if the value is invalid or unsupported.
     */
    DBNode::Value parseValue(const std::string& input, size_t& index) const;

    /**
     * @brief Parses a JSON number (integer or floating point).
     *
     * Accepts the typical JSON number syntax:
     *   - optional minus sign,
     *   - digits,
     *   - optional decimal part,
     *   - optional exponent (e/E with optional sign).
     *
     * @param input JSON string.
     * @param index Reference to the current position; advanced past
     *              the end of the numeric token.
     * @return Parsed double value.
     *
     * @throws std::runtime_error if the token cannot be converted to double.
     */
    double parseNumber(const std::string& input, size_t& index) const;

    /**
     * @brief Skips whitespace characters in the input string.
     *
     * Advances @p index while input[index] is a whitespace character
     * as reported by std::isspace.
     *
     * @param input JSON string.
     * @param index Reference to the current position; updated in place.
     */
    void skipWhitespace(const std::string& input, size_t& index) const;

    /**
     * @brief Parses a JSON array starting at the current position.
     *
     * Expects input[index] to point to '['. Elements are parsed using
     * parseValue() and stored in a DBNode with numeric string keys:
     *   0, 1, 2, ...
     *
     * On return, @p index is advanced past the closing ']'.
     *
     * @param input JSON string.
     * @param index Reference to the current position; advanced as parsing proceeds.
     * @return Shared pointer to the DBNode representing the array.
     *
     * @throws std::runtime_error on invalid syntax.
     */
    std::shared_ptr<DBNode> parseArray(const std::string& input, size_t& index) const;

};

#endif // PARSER_H