#ifndef __YAMLPARSER_H__
#define __YAMLPARSER_H__

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
 * @brief YAML implementation of the IParser class.
 */
class YAMLParser : public IParser {
public:
    std::shared_ptr<Node> parse(const std::string& input) const override;
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override;
    std::shared_ptr<Node> readFromFile(const std::string& filename) const override;

private:
    size_t countLeadingSpaces(const std::string& line) const;
    bool parseLine(const std::string& line, std::string& key, std::string& value) const;
    void adjustIndentation(size_t indent, std::vector<std::shared_ptr<Node>>& nodeStack, std::vector<int>& indentStack) const;
    void trim(std::string& str) const;
    bool isListNode(const std::shared_ptr<Node>& node) const;
    double parseNumber(const std::string& value) const;
    Node::Value parseValue(const std::string& value) const;
};

#endif // __YAMLPARSER_H__
