#ifndef __IPARSER_H__
#define __IPARSER_H__

#include <memory>
#include "DBNode.h"

/**
 * @brief Abstract base class for data parsers.
 */
class IParser {
    public:
        virtual ~IParser() = default;
    
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

#endif // __IPARSER_H__
