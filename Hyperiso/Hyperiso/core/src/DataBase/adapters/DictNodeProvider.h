#ifndef DICTNODEPROVIDER_H
#define DICTNODEPROVIDER_H

#include "INodeProvider.h"
#include "DBManager.h"

/**
 * @file DictNodeProvider.h
 * @brief INodeProvider implementation for dictionary-like formats (JSON/YAML).
 *
 * DictNodeProvider uses DBManager to parse JSON or YAML files into a Node tree.
 */

/**
 * @class DictNodeProvider
 * @brief Node provider for JSON/YAML configuration files.
 *
 * This provider is the counterpart of LhaNodeProvider for dictionary-based
 * formats such as JSON and YAML.
 */
class DictNodeProvider : public INodeProvider {
public:
    /**
     * @brief Constructs a provider bound to a specific JSON/YAML file path.
     *
     * @param src_path Path to the input configuration file.
     */
    DictNodeProvider(fs::path src_path_) : INodeProvider(src_path_) {}

    /**
     * @brief Loads the underlying file and returns it as a Node tree.
     *
     * Internally this delegates to DBManager::read_from_file(), which in turn
     * selects the appropriate parser based on the file extension.
     *
     * @return Shared pointer to the root Node representing the file content.
     */
    std::shared_ptr<Node> provide_db_as_node() override;
};

#endif // DICTNODEPROVIDER_H
