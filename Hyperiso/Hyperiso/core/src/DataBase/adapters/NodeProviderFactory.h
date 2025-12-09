#ifndef NODEPROVIDERFACTORY_H
#define NODEPROVIDERFACTORY_H

#include <memory>
#include "INodeProvider.h"
#include "DictNodeProvider.h"
#include "LhaNodeProvider.h"
#include "DBManager.h"

/**
 * @file NodeProviderFactory.h
 * @brief Factory for creating INodeProvider instances from a file path.
 *
 * This factory inspects the file extension and chooses the appropriate
 * INodeProvider implementation:
 *   - DictNodeProvider for JSON / YAML files,
 *   - LhaNodeProvider  for LHA / SLHA / FLHA files.
 */

/**
 * @class NodeProviderFactory
 * @brief Central factory for building concrete INodeProvider objects.
 *
 * NodeProviderFactory hides the logic of choosing which provider to use
 * for a given source path. It relies on DBManager::EXTENSIONS to map
 * file extensions to parser families.
 */
class NodeProviderFactory {
public:
    /**
     * @brief High-level type of node provider.
     *
     * Currently two families are supported:
     *   - DICT: dictionary-style formats (JSON / YAML),
     *   - LHA : LHA / SLHA / FLHA text files.
     */
    enum class Type { DICT, LHA };

    /**
     * @brief Creates an INodeProvider for a given file path.
     *
     * The extension of @p src_path is inspected and compared against the
     * extension sets registered in DBManager::EXTENSIONS:
     *   - JSON or YAML -> DictNodeProvider,
     *   - LHA / SLHA / FLHA -> LhaNodeProvider.
     *
     * If the extension cannot be recognized, a fatal error is logged.
     *
     * @param src_path Path to the source file whose contents should be exposed
     *                 via an INodeProvider.
     * @return Shared pointer to a concrete INodeProvider.
     */
    static std::shared_ptr<INodeProvider> createNodeProvider(fs::path src_path);
};

#endif // NODEPROVIDERFACTORY_H
