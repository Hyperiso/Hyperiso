#ifndef IFREEZER_H
#define IFREEZER_H

/**
 * @file IFreezer.h
 * @brief Interface for freezing and unfreezing blocks or parameters.
 *
 * This header declares the templated IFreezer interface, which formalizes
 * the contract for freezing / unfreezing:
 *  - whole blocks identified by a (type, name) pair
 *  - individual objects (typically parameters) identified by an ID type.
 *
 * The interface is expressed via deleted static functions so that concrete
 * implementations (e.g. Freezer) provide their own static API without
 * being accidentally instantiated.
 */

/**
 * @example freezer_example.cpp
 * @brief Practical example showing how to freeze and unfreeze parameters and blocks using Freezer.
 *
 * @defgroup FreezeControlModule Freeze/Unfreeze Control System
 * @brief Provides interfaces and classes to freeze and unfreeze parameters and blocks.
 *
 * ## Related Classes
 * - @ref IFreezer
 * - @ref Freezer
 */

/**
 * @class IFreezer
 * @ingroup FreezeControlModule
 * @brief Interface for freezing and unfreezing blocks or parameters.
 *
 * @tparam T Type used to describe the parameter “category” (e.g. ParameterType).
 * @tparam U Type used to identify a block (e.g. std::string).
 * @tparam V Type used to identify a single entity (e.g. ParamId).
 *
 * This interface declares a static API (via deleted functions) intended to be
 * implemented by concrete helper classes such as Freezer. The use of deleted
 * functions prevents accidental usage of the interface directly while still
 * documenting the expected signatures.
 */
template<typename T, typename U, typename V>
class IFreezer {
public:
    /// Virtual destructor for proper polymorphic cleanup if ever inherited non-statically.
    virtual ~IFreezer() = default;

    /**
     * @brief Freezes a block given its type and name.
     *
     * Concrete implementations should prevent further updates to all parameters
     * contained in the block designated by (T, U).
     */
    static void freeze(const T&, const U&) = delete;

    /**
     * @brief Freezes a single object (typically a parameter).
     *
     * Concrete implementations should prevent further updates to the entity
     * identified by V.
     */
    static void freeze(const V&) = delete;

    /**
     * @brief Unfreezes a block given its type and name.
     *
     * Concrete implementations should allow updates again for all parameters
     * contained in the block designated by (T, U).
     */
    static void unfreeze(const T&, const U&) = delete;

    /**
     * @brief Unfreezes a single object (typically a parameter).
     *
     * Concrete implementations should allow updates again for the entity
     * identified by V.
     */
    static void unfreeze(const V&) = delete;
};


#endif // IFREEZER_H
