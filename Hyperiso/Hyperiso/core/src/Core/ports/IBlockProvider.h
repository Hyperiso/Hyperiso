#ifndef IBLOCK_PROVIDER_H
#define IBLOCK_PROVIDER_H

/**
 * @file IBlockProvider.h
 * @brief Defines a generic interface to query and log parameter blocks.
 *
 * This file declares the template interface IBlockProvider, which is meant to be
 * specialized for different block “types” (e.g. ParameterType) and block-name types
 * (e.g. std::string). Concrete implementations (such as BlockProvider) typically
 * delegate to Parameters / BlockAccessor to:
 *
 * - test the existence of a given block,
 * - log the list of available blocks for a given type,
 * - log the content of a specific block.
 */

 /**
 * @class IBlockProvider
 * @brief Generic interface for providing information about parameter blocks.
 *
 * @tparam T Type used to represent a category or model (e.g. ParameterType).
 * @tparam U Type used to represent a block name (e.g. std::string or std::string_view).
 *
 * Implementations are expected to:
 * - check whether a block exists for a given category,
 * - log all available blocks for a category,
 * - log a specific block of a category.
 */
template<typename T, typename U>
class IBlockProvider {
public:
    /**
     * @brief Checks whether a given block exists for a given category/model.
     *
     * @param blockname Name of the block to check.
     * @param type Category / model identifier.
     * @return true if the block exists for this category, false otherwise.
     */
    virtual bool exists(U blockname, T) = 0;

    /**
     * @brief Logs all blocks for a given category/model.
     *
     * Implementations are free to choose how logging is performed
     * (e.g. via Logger, std::cout, etc.).
     *
     * @param type Category / model identifier.
     */
    virtual void log_all_blocks(T type) = 0;

    /**
     * @brief Logs the content of a specific block for a given category/model.
     *
     * @param type Category / model identifier.
     * @param blockname Name of the block to log.
     */
    virtual void log_block(T type, U blockname) = 0;

    /**
     * @brief Retrieve the content of a specific block for a given category/model.
     *
     * @param type Category / model identifier.
     * @param blockname Name of the block to log.
     * @return std::map<LhaID, scalar_t>, The content of the block.
     */
    virtual std::map<LhaID, scalar_t> get_block(ParameterType type, const std::string& blockname) = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~IBlockProvider() = default;
};

#endif