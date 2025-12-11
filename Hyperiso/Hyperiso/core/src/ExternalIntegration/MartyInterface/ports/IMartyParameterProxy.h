#ifndef IMARTYPARAMETERPROXY_H
#define IMARTYPARAMETERPROXY_H

#include "scalar.h"

/**
 * @file IMartyParameterProxy.h
 * @brief Declares a generic interface to access parameters for MARTY adapters.
 *
 * This abstraction allows code-generation utilities to retrieve parameter
 * values without directly depending on the underlying parameter storage
 * (e.g. Hyperiso’s ::Parameters system).
 */

/**
 * @class IMartyParameterProxy
 * @ingroup CodeGenerationModule
 * @brief Interface for accessing parameter values from different backends.
 *
 * @tparam T Type used to identify a block (e.g. ::std::string).
 * @tparam V Type used to identify an entry inside a block (e.g. ::LhaID).
 *
 * Implementations (such as ::MartyParameterProxy) must provide an
 * `operator()(T, V)` returning a scalar value.
 */
template<typename T, typename V>
class IMartyParameterProxy {
public:

    /**
     * @brief Retrieves the numerical value of a parameter.
     *
     * @param x Block identifier.
     * @param y Entry identifier inside the block (e.g. an LHA ID).
     * @return The parameter value as a scalar.
     */
    virtual scalar_t operator()(const T& x, const V& y) const = 0 ;
};

#endif // IMARTYPARAMETERPROXY_H
