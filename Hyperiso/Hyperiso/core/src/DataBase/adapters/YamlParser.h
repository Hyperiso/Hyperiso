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
    bool parseLine(const std::string& line, std::string& key, std::string& value) const;
    void adjustIndentation(size_t indent, std::vector<std::shared_ptr<Node>>& nodeStack, std::vector<int>& indentStack) const;
    static std::shared_ptr<Node> parseYAMLNode(std::istringstream& stream, int indentLevel);
    static void processKeyValue(const std::string& line, std::shared_ptr<Node>& node, std::istringstream& stream, int indentLevel);
    static void processList(std::shared_ptr<Node>& node, std::istringstream& stream, int indentLevel, std::string firstLine = "");
    static std::string trim(const std::string& str);
    static int countLeadingSpaces(const std::string& str);
    bool isListNode(const std::shared_ptr<Node>& node) const;
    double parseNumber(const std::string& value) const;
    static Node::Value parseScalar(const std::string& value);
    Node::Value parseValue(const std::string& value) const;
    static bool isDouble(const std::string& str);
    static bool isInteger(const std::string& str);
};

#endif // __YAMLPARSER_H__
