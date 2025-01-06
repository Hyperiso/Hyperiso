#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <variant>
#include <initializer_list>
#include <stdexcept>

/**
 * @brief A hierarchical data structure with JSON and YAML serialization capabilities.
 */
class Node {
public:
    using Value = std::variant<std::string, int, double, bool, std::shared_ptr<Node>>;

    /**
     * @brief Default constructor for Node.
     */
    Node();

    /**
     * @brief Retrieves a value from the node using a series of keys.
     * @tparam Keys The types of the keys.
     * @param keys The keys specifying the path to the value.
     * @return The retrieved value.
     * @throws std::runtime_error If the key path does not exist.
     */
    template <typename... Keys>
    Value get(Keys&&... keys) const;

    /**
     * @brief Sets a value in the node using a series of keys.
     * @tparam T The type of the value.
     * @tparam Key The type of the first key.
     * @tparam Rest The types of additional keys.
     * @param value The value to set.
     * @param key The first key.
     * @param rest Additional keys specifying the path.
     */
    template <typename T, typename Key, typename... Rest>
    void set(T value, Key&& key, Rest&&... rest);

    /**
     * @brief Retrieves a group of values as a map using a vector of keys.
     * @param keys The keys specifying the path to the group.
     * @return A map of key-value pairs representing the group.
     * @throws std::runtime_error If the key path does not exist.
     */
    std::map<std::string, Value> getGroup(const std::vector<std::string>& keys) const;

    /**
     * @brief Sets a group of values using a vector of keys.
     * @param keys The keys specifying the path to the group.
     * @param groupData A map of key-value pairs representing the group to set.
     */
    void setGroup(const std::vector<std::string>& keys, const std::map<std::string, Value>& groupData);

    /**
     * @brief Prints the node in JSON format to standard output.
     * @param level The indentation level.
     */
    void printJSON(int level = 0) const;

    /**
     * @brief Prints the node in JSON format to a given output stream.
     * @param os The output stream.
     * @param level The indentation level.
     */
    void printJSONToStream(std::ostream& os, int level = 0) const;

    /**
     * @brief Prints the node in YAML format to standard output.
     * @param level The indentation level.
     */
    void printYAML(int level = 0) const;

private:
    std::map<std::string, Value> data_;

    template <typename Key, typename... Rest>
    static Value getRecursive(const std::map<std::string, Value>& map, Key&& key, Rest&&... rest);

    void printValue(const Value& value, int level) const;
    void printValueToStream(std::ostream& os, const Value& value, int level) const;
    void printScalarYAML(const Value& value) const;
};

#ifndef PARSER_H
#define PARSER_H

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
};

/**
 * @brief YAML implementation of the Parser class.
 */
class YAMLParser : public Parser {
public:
    std::shared_ptr<Node> parse(const std::string& input) const override;
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const override;
    std::shared_ptr<Node> readFromFile(const std::string& filename) const override;
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

#endif // NODE_H
