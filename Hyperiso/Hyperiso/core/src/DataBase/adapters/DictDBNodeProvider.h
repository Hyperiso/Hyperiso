#ifndef DICTNODEPROVIDER_H
#define DICTNODEPROVIDER_H

#include "IDBNodeProvider.h"
#include "DBManager.h"

/**
 * @file DictDBNodeProvider.h
 * @brief IDBNodeProvider implementation for dictionary-like formats (JSON/YAML).
 *
 * DictDBNodeProvider uses DBManager to parse JSON or YAML files into a DBNode tree.
 */

/**
 * @class DictDBNodeProvider
 * @brief DBNode provider for JSON/YAML configuration files.
 *
 * This provider is the counterpart of LhaDBNodeProvider for dictionary-based
 * formats such as JSON and YAML.
 */
class DictDBNodeProvider : public IDBNodeProvider {
public:
    /**
     * @brief Constructs a provider bound to a specific JSON/YAML file path.
     *
     * @param src_path Path to the input configuration file.
     */
    DictDBNodeProvider(fs::path src_path_) : IDBNodeProvider(src_path_) {}

    /**
     * @brief Loads the underlying file and returns it as a DBNode tree.
     *
     * Internally this delegates to DBManager::read_from_file(), which in turn
     * selects the appropriate parser based on the file extension.
     *
     * @return Shared pointer to the root DBNode representing the file content.
     */
    std::shared_ptr<DBNode> provide_db_as_node() override;
};

#endif // DICTNODEPROVIDER_H
