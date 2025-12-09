#ifndef IDATALOADER_H
#define IDATALOADER_H

#include <memory>
#include <filesystem>

/**
 * @file IDataLoader.h
 * @brief Defines a generic interface for data loading operations.
 *
 * This file declares the abstract template class IDataLoader, which provides
 * a common interface for loading data from a file into a destination object.
 */

/**
 * @defgroup DataLoadersModule Data Loading System
 * @brief System for loading parameters and correlations from external files.
 *
 * This module defines classes responsible for reading structured files (JSON, YAML, etc.)
 * and filling memory structures like:
 * - BlockAccessor (for parameters)
 * - CorrelationMatrixPair (for correlations)
 *
 * ## Overview
 * The DataLoader system is designed to provide a clean and flexible way to load data into memory:
 *
 * - **IDataLoader** is a generic interface for loaders.
 * - **ParamBlockLoader** loads parameter blocks into BlockAccessor.
 * - **CorrelationLoader** loads statistical and systematic correlations into CorrelationMatrixPair.
 *
 * ## Related Classes
 * - @ref IDataLoader
 * - @ref ParamBlockLoader
 * - @ref CorrelationLoader
 * ## Diagram
 * @dot
 * digraph DataLoaders {
 *   node [shape=record, fontname=Helvetica, fontsize=10];
 *
 *   IDataLoader [label="{ IDataLoader<T> | load() }"];
 *   ParamBlockLoader [label="{ ParamBlockLoader | load() }"];
 *   CorrelationLoader [label="{ CorrelationLoader<T> | load() }"];
 *
 *   IDataLoader -> ParamBlockLoader;
 *   IDataLoader -> CorrelationLoader;
 * }
 * @enddot
 */

namespace fs = std::filesystem;

/**
 * @class IDataLoader
 * @ingroup DataLoadersModule
 * @brief Abstract interface for loading data into an object.
 *
 * @tparam T Type of object into which the data will be loaded.
 */
template<typename T>
class IDataLoader {
public:

    /**
     * @brief Loads data from a file into a destination object.
     *
     * @param dest     Shared pointer to the destination object.
     * @param src_file Path to the source file to load from.
     */
    virtual void load(std::shared_ptr<T> dest, fs::path src_file) = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~IDataLoader() = default;
};

#endif // IDATALOADER_H
