#ifndef DBNODE_H
#define DBNODE_H

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
     * @brief Retrieves all the stored keys below the current node.
     * @return A vector of all the stored keys below the current node.
     */
    std::vector<std::string> get_keys();

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

#endif