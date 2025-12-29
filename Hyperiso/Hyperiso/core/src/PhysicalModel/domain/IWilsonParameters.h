#ifndef IWILSON_PARAMETERS_H
#define IWILSON_PARAMETERS_H

#include <memory>

#include "wgroup_ids.hpp"

class IBlockComposer;

/**
 * @file IWilsonParameterHelper.h
 * @brief Interface for creating and managing Wilson “helper” parameter blocks.
 *
 * Wilson computations rely on a set of auxiliary blocks:
 *  - scale-independent parameters (generation-dependent constants, etc.),
 *  - matching-scale parameters (mu_W dependent quantities),
 *  - running-scale parameters (eta factors, evolution matrices, etc.).
 *
 * Implementations of this interface are responsible for composing those blocks
 * via an injected @ref IBlockComposer (dependency graph / dependent blocks system),
 * and for cleaning up if needed.
 *
 * Typical usage:
 * @code
 *   auto helper = std::make_shared<WilsonParameterHelper>(iblock_c);
 *   helper->init(gen, GroupMapper::to_id(WGroup::B));
 *   ...
 *   helper->cleanup();
 * @endcode
 *
 * @see IBlockComposer
 * @see DependentBlock
 * @see WilsonParameterHelper
 */
class IWilsonParameterHelper {
public:
    /**
     * @brief Constructs the helper using a block composer.
     * @param iblock_c Block composer used to register dependent blocks/params.
     */
    IWilsonParameterHelper(std::shared_ptr<IBlockComposer> iblock_c) : iblock_c(iblock_c) {}

    virtual ~IWilsonParameterHelper() = default;

    /**
     * @brief Initializes all required helper blocks.
     *
     * Implementations should be idempotent: calling init() twice should not
     * register duplicates.
     *
     * @param gen Generation index (used by some blocks, e.g. lepton mass choice).
     * @param grp Wilson group identifier (may steer which running matrices are built).
     */
    virtual void init(int gen, WGroupId grp) = 0;

    /**
     * @brief Cleans up internal state and/or unregisters dependent blocks (if applicable).
     *
     * Some implementations only reset internal flags, others may remove composed blocks.
     */
    virtual void cleanup() = 0;

    /**
     * @brief Returns whether init() has already been successfully called.
     */
    bool is_init() const {return initialized;}

protected:
    /**
     * @brief Builds scale-independent Wilson helper blocks (generation-dependent).
     */
    virtual void init_scale_independent_block(int gen) = 0;

    /**
     * @brief Builds matching-scale dependent Wilson helper blocks.
     */
    virtual void init_matching_block() = 0;

    /**
     * @brief Builds running-scale dependent Wilson helper blocks.
     */
    virtual void init_running_block(WGroupId grp) = 0;

    /// Block composer (dependency engine) used to register dependent blocks.
    std::shared_ptr<IBlockComposer> iblock_c;

    /// Tracks whether the helper has already been initialized.
    bool initialized{false};
};

#endif