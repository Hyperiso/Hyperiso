#ifndef IMARTY_WILSON_PROXY_H
#define IMARTY_WILSON_PROXY_H

#include <set>
#include <unordered_set>
#include <string>

/**
 * @file IMartyWilsonProxy.h
 * @brief Interface for MARTY-based Wilson coefficient calculations.
 *
 * This header defines @ref IMartyWilsonProxy, an abstract interface
 * exposing high-level operations related to Wilson coefficient
 * generation and dependency tracking using MARTY.
 *
 * @tparam T Type used to represent interpreted parameter dependencies
 *           (typically @ref InterpretedParam).
 */

/**
 * @class IMartyWilsonProxy
 * @ingroup MartyIntegrationModule
 * @brief Abstract proxy for Wilson coefficient calculations via MARTY.
 *
 * This interface defines the contract for objects responsible for:
 *  - generating Wilson coefficients,
 *  - managing numerical libraries,
 *  - tracking parameter dependencies,
 *  - exposing special blocks required for calculations.
 *
 * Implementations typically wrap a concrete adapter such as
 * @ref MartyWilsonAdapter.
 *
 * @tparam T Type representing a resolved/interpreted parameter.
 */
template <typename T>
class IMartyWilsonProxy {
public:
    /**
     * @brief Runs the full Wilson coefficient calculation pipeline.
     *
     * @param wilson      Name of the Wilson coefficient.
     * @param model       Model name (SM, THDM, MSSM, etc.).
     * @param Q_match     Matching scale.
     * @param model_path  Path to the model definition file.
     * @param new_params  Whether parameters should be regenerated.
     */
    virtual void calculate(std::string wilson, std::string model, double Q_match, std::string model_path, bool new_params = false) = 0;

    /**
     * @brief Returns the set of special parameter blocks.
     *
     * These blocks require custom treatment during parameter extraction
     * and numerical code generation.
     *
     * @return Set of block names.
     */
    virtual std::set<std::string>  get_special_blocks() = 0;

    /**
     * @brief Returns the set of parameter dependencies for a given Wilson coefficient.
     *
     * @param wilson Name of the Wilson coefficient.
     * @return Set of interpreted parameters required by this Wilson coefficient.
     */
    virtual std::unordered_set<T> get_dependencies(std::string wilson) = 0;

};

#endif // IMARTY_WILSON_PROXY_H