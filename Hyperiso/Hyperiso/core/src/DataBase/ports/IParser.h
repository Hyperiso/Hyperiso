#ifndef IPARSER_H
#define IPARSER_H

#include <memory>
#include "DBNode.h"

/**
 * @file IParser.h
 * @brief Abstract interface for (de)serializing data as Node trees.
 *
 * An IParser provides a common API to:
 *   - parse a textual representation into a Node hierarchy,
 *   - serialize a Node hierarchy to a file,
 *   - read a Node hierarchy back from a file.
 *
 * Concrete implementations (e.g. YAMLParser, JSONParser, LhaParser)
 * specify the actual format.
 */

/**
 * @class IParser
 * @brief Abstract base class for data parsers.
 *
 * This interface is format-agnostic. It exposes three core operations:
 *   - parse(): build a Node tree from an in-memory string,
 *   - writeToFile(): serialize a Node tree to a file,
 *   - readFromFile(): read a Node tree from a file (and typically
 *     delegate to parse()).
 */
class IParser {
    public:
        /// Virtual destructor to allow proper cleanup via base pointer.
        virtual ~IParser() = default;
    
        /**
         * @brief Parses an input string into a Node hierarchy.
         *
         * The exact syntax and semantics of the input depend on the
         * concrete implementation (YAML, JSON, LHA, ...).
         *
         * @param input Text buffer to parse.
         * @return Shared pointer to the root Node of the parsed tree.
         *
         * @throws std::runtime_error (or derived) on parsing errors.
         */
        virtual std::shared_ptr<Node> parse(const std::string& input) const = 0;
    
        /**
         * @brief Serializes a Node hierarchy to a file.
         *
         * Implementations are free to choose their own textual encoding
         * (YAML, JSON, etc.) but should be compatible with readFromFile()
         * for round-trip usage.
         *
         * @param filename Path of the target file.
         * @param root     Root Node to serialize.
         *
         * @throws std::runtime_error if the file cannot be opened or a
         *         write error occurs.
         */
        virtual void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const = 0;
    
        /**
         * @brief Reads and parses a Node hierarchy from a file.
         *
         * Typical implementations:
         *   - read the entire file into a string,
         *   - call parse() on that string.
         *
         * @param filename Path of the file to read.
         * @return Shared pointer to the root Node of the parsed tree.
         *
         * @throws std::runtime_error if the file cannot be opened or if
         *         parsing fails.
         */
        virtual std::shared_ptr<Node> readFromFile(const std::string& filename) const = 0;
    };

#endif // IPARSER_H
