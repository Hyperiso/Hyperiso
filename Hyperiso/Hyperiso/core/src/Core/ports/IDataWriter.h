#ifndef IDATAWRITER_H
#define IDATAWRITER_H

#include <memory>
#include <filesystem>

/**
 * @file IDataWriter.h
 * @brief Defines a generic interface for data writing operations.
 *
 * This file declares the abstract template class IDataWriter, which provides
 * a common interface for writing data from a BlockAccessor into a destination object.
 */

/**
 * @defgroup DataWritersModule Data Writing System
 * @brief System for writing parameters and correlations to external files (along with the DataBase layer).
 *
 * This module defines classes responsible for writing BlockAccessor to structured files (JSON, YAML, etc.)
 *
 * ## Overview
 * The DataWriter system is designed to provide a clean and flexible way to write data into memory:
 *
 * - **IDataWriter** is a generic interface for writers.
 * - **ParamBlockWriter** Write parameter blocks from a BlockAccessor to a DBNode.
 *
 * ## Related Classes
 * - @ref IDataWriter
 * - @ref ParamBlockWriter
 * ## Diagram
 * @dot
 * digraph DataWriters {
 *   node [shape=record, fontname=Helvetica, fontsize=10];
 *
 *   IDataWriter [label="{ IDataLoader<T> | write() }"];
 *   ParamBlockWriter [label="{ ParamBlockLoader | write() }"];
 *
 *   IDataWriter -> ParamBlockWriter;
 * }
 * @enddot
 */

namespace fs = std::filesystem;

/**
 * @class IDataWriter
 * @ingroup DataLoadersModule
 * @brief Abstract interface for writing data into an object.
 *
 * @tparam T Type of object into which the data will be written.
 *  * @tparam P Type of object from which the data will be written.
 */
template<typename T, typename P>
class IDataWriter {
public:

    /**
     * @brief Write data from BlockAccessor into a DBNode.
     *
     * @param dest     the destination object.
     * @param src      src object (BlockAccessor).
     */
    virtual void write(T dest, P src) = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~IDataWriter() = default;
};

#endif // IDATAWRITER_H
