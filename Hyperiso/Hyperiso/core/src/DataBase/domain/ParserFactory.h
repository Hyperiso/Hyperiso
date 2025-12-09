#ifndef PARSERFACTORY_H
#define PARSERFACTORY_H

#include <memory>
#include "IParser.h"
#include "JsonParser.h"
#include "YamlParser.h"
#include "LhaParser.h"

/**
 * @brief Factory for creating IParser instances.
 */
class ParserFactory {
public:
    enum class Type { JSON, YAML, LHA };

    /**
     * @brief Creates a parser instance of the specified type.
     * @param type The type of parser to create.
     * @return A unique pointer to the created parser.
     */
    static std::shared_ptr<IParser> createParser(Type type);
};

inline std::shared_ptr<IParser> ParserFactory::createParser(Type type) {
    switch (type) {
    case Type::JSON:
        return std::make_shared<JSONParser>();
    case Type::YAML:
        return std::make_shared<YAMLParser>();
    case Type::LHA:
        return std::make_shared<LhaParser>();
    default:
        throw std::invalid_argument("Unknown parser type");
    }
}

#endif // PARSERFACTORY_H