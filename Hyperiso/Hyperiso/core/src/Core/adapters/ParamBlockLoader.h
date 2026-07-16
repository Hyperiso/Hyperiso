#ifndef BLOCKSCREATOR_H
#define BLOCKSCREATOR_H

#include <memory>

#include "IDataLoader.h"
#include "ILhaPrototypeRegistry.h"
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
class ParamBlockLoader : public IDataLoader<BlockAccessor>, public ILhaPrototypeRegistry {
public:
    /**
     * @brief Loads parameter blocks from a file into the given BlockAccessor.
     *
     * @param dest     Shared pointer to destination BlockAccessor.
     * @param src_file Path to the source file.
     */
    void load(std::shared_ptr<BlockAccessor> dest, fs::path src_file, bool block_in_blocks=false) override;

    /**
     * @brief Registers an additional LHA prototype to apply before parsing LHA files.
     */
    void add_lha_prototype(BlockName blockName,
                           size_t itemCount=2,
                           size_t valueIdx=1,
                           int scaleIdx=-1,
                           int rgIdx=-1,
                           int binIdx=-1,
                           bool globalScale=false) override;

    /**
     * @brief Registers several additional LHA prototypes.
     */
    void add_lha_prototypes(const std::vector<LhaPrototypeSpec>& prototypes) override;

private:
    void apply_lha_prototypes(std::shared_ptr<IDBNodeProvider> provider) const;

    std::vector<LhaPrototypeSpec> lha_prototypes;
};

#endif // BLOCKSCREATOR_H
