#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <memory>

#include "Include.h"
#include "IDataWriter.h"
#include "BlockAccessor.h"
#include "DBNode.h"
#include "MemoryManager.h"
#include "Parameters.h"
#include "ParamBlockWriter.h"
#include "JsonParser.h"
#include "YamlParser.h"
#include "LhaParser.h"
#include "Include.h"

/**
 * @file FileWriter.h
 * @brief Write parameter blocks from a BlockAccessor into a file (LHA, JSON and YAML format).
 *
 * This file defines the FileWriter class, which specializes the IDataWriter interface
 * to write parameter blocks into a file with differents formats (.json, .yaml, .lha, .slha and .flha).
 */

/**
 * @class FileWriter
 * @ingroup DataWritersModule
 * @brief Writer class to populate a file from a BlockAccessor.
 */
class FileWriter : public IDataWriter<const std::string&, std::shared_ptr<BlockAccessor>> {
public:
    /**
     * @brief Write node parameters from a BlockAccessor into the destination file. If no BlockAccessor is provided, use the BlockAccessor from memory manager.
     *
     * @param dest     path to the destination file.
     * @param mode     Unused parameter.
     */
    void write(const std::string& dest, std::shared_ptr<BlockAccessor> src = nullptr) override;
};

#endif // FILE_WRITER_H
