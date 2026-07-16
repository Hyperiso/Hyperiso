#ifndef ABSTRACTCONFIG_H
#define ABSTRACTCONFIG_H

/**
 * @file AbstractConfig.h
 * @brief Base type for lightweight configuration structures.
 *
 * This header defines the AbstractConfig base struct, from which all
 * configuration POD-like types (such as WilsonBuildConfig, WilsonRequest,
 * AlphasConfig, MassConfig, etc.) inherit.
 *
 * The main purpose of this type is to provide a common polymorphic base
 * with a virtual destructor, ensuring that configuration objects can be
 * handled via base-class pointers and safely destroyed.
 */

/**
 * @struct AbstractConfig
 * @brief Polymorphic base class for configuration types.
 *
 * AbstractConfig is a minimal base type that introduces a virtual
 * destructor, allowing derived configuration structures to be managed
 * through pointers or references to AbstractConfig without risking
 * undefined behaviour at destruction time.
 *
 * It does not impose any interface beyond lifetime management and is
 * intended to keep configuration types lightweight and POD-like
 * wherever possible.
 */
struct AbstractConfig {
    /**
     * @brief Virtual destructor to allow safe polymorphic deletion.
     *
     * Declaring a virtual destructor ensures that when a configuration
     * object is deleted through a pointer to AbstractConfig, the
     * destructor of the most-derived type is invoked.
     */
    virtual ~AbstractConfig() = default;
};

#endif // ABSTRACTCONFIG_H
