#ifndef INODEPROVIDER_H
#define INODEPROVIDER_H

#include <memory>
#include <filesystem>

#include "DBNode.h"

namespace fs = std::filesystem;

/**
 * @file INodeProvider.h
 * @brief Interface for objects that expose a database/configuration as a Node tree.
 *
 * Node providers are small adaptors that know how to obtain a Node
 * representation of some backing data source (file, DB, etc.).
 */

/**
 * @class INodeProvider
 * @brief Abstract provider of Node-based views of external data.
 *
 * An INodeProvider typically:
 *   - stores the location of the underlying data (e.g. a file path),
 *   - offers a single method to materialize that data as a Node tree.
 *
 * Concrete implementations in this project:
 *   - DictNodeProvider: for JSON/YAML configuration files,
 *   - LhaNodeProvider: for LHA/SLHA/FLHA input files.
 */
class INodeProvider {
protected:
    fs::path src_path;

public:
    /**
     * @brief Constructs a provider bound to a given source path.
     *
     * @param src_path Location of the backing resource to read from.
     */
    INodeProvider(fs::path src_path) : src_path(src_path) {}
    
    virtual ~INodeProvider() = default;

    /**
     * @brief Returns the underlying data as a Node tree.
     *
     * Implementations are free to choose how they obtain the Node:
     *   - parsing a file (via DBManager),
     *   - querying a database,
     *   - composing multiple sources, etc.
     *
     * @return Shared pointer to the root Node representing the data.
     *
     * @throws std::runtime_error or implementation-specific exceptions
     *         in case of I/O or parsing failures.
     */
    virtual std::shared_ptr<Node> provide_db_as_node() = 0; 
};

#endif // INODEPROVIDER_H
