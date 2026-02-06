#ifndef BLOCKSCREATOR_H
#define BLOCKSCREATOR_H

#include <memory>

#include "IDataLoader.h"
#include "BlockAccessor.h"
#include "DBNodeProviderFactory.h"

/**
 * @file ParamBlockLoader.h
 * @brief Loads parameter blocks from a file into a BlockAccessor.
 *
 * This file defines the ParamBlockLoader class, which specializes the IDataLoader interface
 * to load parameter blocks into a BlockAccessor instance.
 */

/**
 * @class ParamBlockLoader
 * @ingroup DataLoadersModule
 * @brief Loader class to populate a BlockAccessor from a source file.
 */
class ParamBlockLoader : public IDataLoader<BlockAccessor> {
public:
    /**
     * @brief Loads parameter blocks from a file into the given BlockAccessor.
     *
     * @param dest     Shared pointer to destination BlockAccessor.
     * @param src_file Path to the source file.
     */
    void load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) override;
};

#endif // BLOCKSCREATOR_H
