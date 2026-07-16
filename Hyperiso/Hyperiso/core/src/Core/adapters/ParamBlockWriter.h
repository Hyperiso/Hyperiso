#ifndef BLOCKS_WRITER_H
#define BLOCKS_WRITER_H

#include <memory>

#include "IDataWriter.h"
#include "BlockAccessor.h"
#include "DBNodeProviderFactory.h"
#include "DBNode.h"

/**
 * @file ParamBlockWriter.h
 * @brief Write parameter blocks from a BlockAccessor into a DBNode.
 *
 * This file defines the ParamBlockWriter class, which specializes the IDataWriter interface
 * to write parameter blocks into a DBNode instance.
 */

/**
 * @class ParamBlockWriter
 * @ingroup DataWritersModule
 * @brief Writer class to populate a DbNode from a BlockAccessor.
 */
class ParamBlockWriter : public IDataWriter<std::shared_ptr<DBNode>, std::shared_ptr<BlockAccessor>> {
public:
    /**
     * @brief Write node parameters from a BlockAccessor into the given destination node.
     *
     * @param dest     Shared pointer to destination BlockAccessor.
     * @param src      Shared pointer to src BlockAccessor.
     */
    void write(std::shared_ptr<DBNode> dest, std::shared_ptr<BlockAccessor> src) override;
};

#endif // BLOCKS_WRITER_H
