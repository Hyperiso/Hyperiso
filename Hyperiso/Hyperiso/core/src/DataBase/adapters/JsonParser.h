#ifndef PARSER_H
#define PARSER_H

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
 * @brief JSON implementation of the IParser class.
 */
class JSONParser : public IParser {
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

#endif // PARSER_H