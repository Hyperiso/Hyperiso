#ifndef LHANODEPROVIDER_H
#define LHANODEPROVIDER_H

#include <string>

#include "IDBNodeProvider.h"
#include "DBManager.h"

/**
 * @file LhaDBNodeProvider.h
 * @brief IDBNodeProvider implementation backed by an LHA/SLHA/FLHA file.
 *
 * LhaDBNodeProvider uses DBManager under the hood to:
 *   - parse the given LHA-like file,
 *   - obtain a DBNode tree representation,
 *   - optionally extend the set of known LHA block prototypes.
 */

/**
 * @class LhaDBNodeProvider
 * @brief DBNode provider for LHA / SLHA / FLHA sources.
 *
 * This provider delegates the heavy lifting to DBManager:
 *   - DBManager::read_from_file() for parsing,
 *   - DBManager::add_lha_prototype() for enriching the set of block
 *     prototypes used by the LHA parser.
 */
class LhaDBNodeProvider : public IDBNodeProvider {
public:
    /**
     * @brief Constructs a provider bound to a specific source path.
     *
     * @param src_path Path to an LHA / SLHA / FLHA file.
     */
    LhaDBNodeProvider(fs::path src_path_) : IDBNodeProvider(src_path_) {}

    /**
     * @brief Loads the underlying LHA file and returns it as a DBNode tree.
     *
     * Internally this is equivalent to:
     * @code
     *   DBManager db;
     *   return db.read_from_file(src_path);
     * @endcode
     *
     * @return Shared pointer to the root DBNode representing the LHA content.
     *
     * @throws std::runtime_error if file or parsing errors occur.
     */
    std::shared_ptr<DBNode> provide_db_as_node() override;

    /**
     * @brief Adds a new LHA block prototype to the global set.
     *
     * This is a convenience wrapper around DBManager::add_lha_prototype(),
     * allowing clients to register additional block layouts before parsing.
     *
     * The arguments follow the Prototype layout:
     *   - @p blockName   : block identifier (case-insensitive),
     *   - @p itemCount   : number of columns in each line,
     *   - @p valueIdx    : index of the value column,
     *   - @p scaleIdx    : index of the scale column (-1 if none),
     *   - @p rgIdx       : index of the renormalization scheme column (-1 if none),
     *   - @p binIdx      : index of the bin lower edge (-1 if unbinned),
     *   - @p globalScale : true if the block uses a global Q= header.
     *
     * Input indices are checked for basic consistency and will trigger
     * a logged error if the configuration is invalid.
     */
    void add_lha_prototype(BlockName blockName, size_t itemCount=2, size_t valueIdx=1, int scaleIdx=-1, int rgIdx=-1, int binIdx=-1, bool globalScale=false);

};

#endif // LHANODEPROVIDER_H
