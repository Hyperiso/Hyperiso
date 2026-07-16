#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <memory>
#include <string>
#include <vector>

#include "BlockAccessor.h"
#include "IDataWriter.h"
#include "ParamID.h"

/**
 * @file FileWriter.h
 * @brief Export the current HyperIso parameter database to JSON, YAML or LHA.
 *
 * The writer serializes parameters as block/id entries. The output format is
 * selected from the destination suffix:
 * - `.json` for JSON,
 * - `.yaml` or `.yml` for YAML,
 * - `.lha`, `.slha` or `.flha` for Les Houches formats.
 *
 * When no explicit BlockAccessor is supplied, the writer exports the current
 * database managed by MemoryManager. Callers may export the complete database,
 * a selected set of blocks, or a selected set of ParamId entries.
 */
class FileWriter : public IDataWriter<const std::string&, std::shared_ptr<BlockAccessor>> {
public:
    /**
     * @brief Export all blocks from an accessor or from the current database.
     *
     * @param dest Destination filename. Its suffix selects the output format.
     * @param src Optional source accessor. When null, the current database is used.
     *
     * @throws std::logic_error If the current database is requested before initialization.
     * @throws std::invalid_argument If the destination suffix is unsupported.
     * @throws std::runtime_error If the destination cannot be written.
     */
    void write(const std::string& dest, std::shared_ptr<BlockAccessor> src = nullptr) override;

    /**
     * @brief Export only the requested blocks.
     *
     * @param dest Destination filename.
     * @param block_names Block names or aliases to export.
     * @param src Optional source accessor. When null, the current database is used.
     *
     * @throws std::invalid_argument If a requested block does not exist.
     */
    void write_blocks(const std::string& dest,
                      const std::vector<BlockName>& block_names,
                      std::shared_ptr<BlockAccessor> src = nullptr);

    /**
     * @brief Export only the requested block/id parameters.
     *
     * ParameterType, when present in a ParamId, is used to resolve the source
     * Parameters repository. Untyped identifiers are resolved from the supplied
     * accessor or from the merged current database.
     *
     * @param dest Destination filename.
     * @param parameter_ids Parameter identifiers to export.
     * @param src Optional source accessor used for untyped identifiers.
     *
     * @throws std::invalid_argument If a requested parameter does not exist.
     */
    void write_parameters(const std::string& dest,
                          const std::vector<ParamId>& parameter_ids,
                          std::shared_ptr<BlockAccessor> src = nullptr);

private:
    static std::shared_ptr<BlockAccessor> resolve_source(std::shared_ptr<BlockAccessor> src);
    static void write_accessor(const std::string& dest, const std::shared_ptr<BlockAccessor>& src);
};

#endif // FILE_WRITER_H
