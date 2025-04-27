#ifndef IFREEZER_H
#define IFREEZER_H

/**
 * @example freezer_example.cpp
 * @brief Practical example showing how to freeze and unfreeze parameters and blocks using Freezer.
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
 * Defines the contract for freezing and unfreezing functionalities.
 */
template<typename T, typename U, typename V>
class IFreezer {
public:
    virtual ~IFreezer() = default;

    /**
     * @brief Freezes a block given its type and name.
     */
    static void freeze(const T&, const U&) = delete;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void freeze(const V&) = delete;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void unfreeze(const T&, const U&) = delete;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void unfreeze(const V&) = delete;
};


#endif // IFREEZER_H
