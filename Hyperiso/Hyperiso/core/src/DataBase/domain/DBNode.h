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
#include <algorithm>

#include "BlockName.h"

/**
 * @file DBNode.h
 * @brief Hierarchical key–value tree with JSON/YAML serialization.
 *
 * Node represents a generic hierarchical data structure:
 *   - keys are BlockName objects,
 *   - values are scalars (int, double, bool, BlockName) or nested Node
 *     instances, including lists of nodes.
 *
 * The class provides:
 *   - typed hierarchical access via get() / set(),
 *   - grouping operations via getGroup() / setGroup(),
 *   - JSON and YAML pretty-printers (to stdout or to an arbitrary stream),
 *   - simple structural queries (contains(), countChildren(), isList()).
 */
class Node {
public:
    /**
     * @brief Variant type used to store values in the tree.
     *
     * A value can be:
     *   - BlockName          : named scalar value,
     *   - int, double, bool  : numeric / boolean scalar,
     *   - std::shared_ptr<Node>                    : nested object,
     *   - std::vector<std::shared_ptr<Node>>       : list of objects
     *                                               (used for arrays).
     */
    using Value = std::variant<BlockName, int, double, bool, std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>>;


    /**
     * @brief Default constructor.
     *
     * Constructs an empty node with no children.
     */
    Node();

    /**
     * @brief Retrieves a value from the node using a sequence of keys.
     *
     * This is a generic hierarchical accessor:
     * @code
     *   Node root;
     *   // ...
     *   auto v = root.get("BLOCK", "SUBBLOCK", "ENTRY");
     * @endcode
     *
     * The provided keys are interpreted as a path: at each step, the
     * corresponding child must be a nested Node (stored in a
     * std::shared_ptr<Node>) until the final element is reached.
     *
     * @tparam Keys Types of the keys (typically BlockName, std::string, const char*).
     * @param keys  Path of keys leading to the desired value.
     * @return The stored Value at the end of the path.
     *
     * @throws std::runtime_error If the path does not exist or the structure
     *         along the way is not compatible (e.g. non-Node where a Node is expected).
     */
    template <typename... Keys>
    Value get(Keys&&... keys) const;

    /**
     * @brief Returns all the keys stored at the current node level.
     *
     * This does not recurse into children; it only returns the immediate
     * keys present in the underlying map.
     *
     * @return Vector of BlockName keys present in this node.
     */
    std::vector<BlockName> get_keys();

    /**
     * @brief Sets a value in the node using a sequence of keys.
     *
     * Intermediate nodes on the path are created as needed and stored as
     * std::shared_ptr<Node>. The last key in the sequence is associated
     * with the provided value.
     *
     * @tparam T    Type of the value to store (compatible with Value).
     * @tparam Key  Type of the first key.
     * @tparam Rest Types of additional keys.
     *
     * @param value Value to store.
     * @param key   First key in the path.
     * @param rest  Remaining keys specifying the path.
     */
    template <typename T, typename Key, typename... Rest>
    void set(T value, Key&& key, Rest&&... rest);

    /**
     * @brief Retrieves a whole sub-tree as a map from a path of keys.
     *
     * The keys vector is interpreted as a path of nested nodes. At the end
     * of the path, the underlying map<BlockName,Value> is returned.
     *
     * @param keys Sequence of keys describing the path to the group node.
     * @return Map of key–value pairs stored at the target node.
     *
     * @throws std::runtime_error If the path does not exist or does not
     *         terminate on a Node.
     */
    std::map<BlockName, Value> getGroup(const std::vector<BlockName>& keys) const;

    /**
     * @brief Sets the content of a whole sub-tree from a map of values.
     *
     * The keys vector is interpreted as a path. Intermediate nodes are
     * created as needed. At the end of the path, the node's internal map
     * is replaced with @p groupData.
     *
     * @param keys      Path to the group node to fill.
     * @param groupData Map of key–value pairs to store at that node.
     */
    void setGroup(const std::vector<BlockName>& keys, const std::map<BlockName, Value>& groupData);

    /**
     * @brief Prints the node in JSON format to std::cout.
     *
     * This produces a pretty-printed JSON-like representation of the
     * current node and all its descendants.
     *
     * @param level Indentation level (used internally for recursion).
     */
    void printJSON(int level = 0) const;

    /**
     * @brief Prints the node in JSON format to a given output stream.
     *
     * @param os    Output stream.
     * @param level Indentation level (used internally for recursion).
     */
    void printJSONToStream(std::ostream& os, int level = 0) const;

    /**
     * @brief Prints the node in YAML format to a given output stream.
     *
     * Lists and nested objects are rendered using standard YAML
     * indentation and list syntax.
     *
     * @param os    Output stream.
     * @param level Indentation level (used internally for recursion).
     */
    void printYAMLToStream(std::ostream& os, int level = 0) const;

    /**
     * @brief Prints the node in YAML format to std::cout.
     *
     * @param level Indentation level (used internally for recursion).
     */
    void printYAML(int level = 0) const;

    /**
     * @brief Checks whether the node contains a direct child with the given key.
     *
     * This does not inspect nested children; it only checks the current
     * level of the map.
     *
     * @param key Key to look up.
     * @return true if @p key exists in this node, false otherwise.
     */
    bool contains(const BlockName& key) const;

    /**
     * @brief Returns the number of direct children of this node.
     *
     * Equivalent to data_.size().
     */
    int countChildren() const;

    /**
     * @brief Returns true if this node represents a "list" node.
     *
     * A node is considered a list node if its first stored value is a
     * std::vector<std::shared_ptr<Node>>. This is a heuristic used
     * internally to differentiate arrays from objects during printing.
     */
    bool isList() const;
private:
    /// Underlying storage for key–value pairs at this node level.
    std::map<BlockName, Value> data_;

    // Recursive lookup helpers used by get(...).

    /**
     * @brief Base case of the recursive getter: single key.
     *
     * Looks up @p key in the given map and returns the associated value.
     * This is the final step in the recursive getRecursive chain.
     */
    template <typename Map, typename Key>
    Value getRecursive(const Map& map, Key&& key) const;

    /**
     * @brief Recursive case of the getter: key followed by more keys.
     *
     * Looks up @p key in @p map, expects the value to be a nested Node
     * (std::shared_ptr<Node>), and continues the lookup on that child
     * with the remaining keys.
     */
    template <typename Map, typename Key, typename Next, typename... Rest>

    /// JSON helper: prints any Value to std::cout with the given indentation.
    Value getRecursive(const Map& map, Key&& key, Next&& next, Rest&&... rest) const;

    /// JSON helper: prints any Value to an arbitrary stream with indentation.
    void printValue(const Value& value, int level) const;

    /// JSON helper: prints any Value to an arbitrary stream with indentation.
    void printValueToStream(std::ostream& os, const Value& value, int level) const;

    // YAML helper: prints a scalar Value to an arbitrary stream.
    void printScalarYAMLToStream(std::ostream& os, const Value& value) const;

    /// YAML helper: prints a scalar Value to std::cout.
    void printScalarYAML(const Value& value) const;

    /**
     * @brief Helper to test if a child node should be treated as a list node.
     *
     * A node is considered a list node if it is non-null, non-empty, and its
     * first stored value is a std::vector<std::shared_ptr<Node>>.
     */
    bool isListNode(const std::shared_ptr<Node>& node) const;
};

#include "DBNode.tpp"

#endif