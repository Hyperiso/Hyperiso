#ifndef IMARTY_WILSON_ADAPTER_H
#define IMARTY_WILSON_ADAPTER_H

#include <string>
#include <set>
#include <unordered_set>

/**
 * @file IMartyWilsonAdapter.h
 * @brief Declares an abstraction layer for MARTY-based Wilson coefficient calculations.
 *
 * The ::IMartyWilsonAdapter interface defines a minimal contract for any
 * adapter that:
 *  - triggers a Wilson-coefficient calculation for a given model,
 *  - exposes special block names,
 *  - exposes discovered parameter dependencies.
 */

/**
 * @class IMartyWilsonAdapter
 * @ingroup CodeGenerationModule
 * @brief Interface for model/Wilson adapters built on MARTY.
 *
 * @tparam T Type used to represent parameter dependencies
 *           (for example ::InterpretedParam).
 *
 * The interface provides three main responsibilities:
 *  - launching a calculation (::calculate),
 *  - querying blocks requiring special treatment (::get_special_blocks),
 *  - exposing the set of discovered dependencies for a given Wilson basis
 *    (::get_dependencies).
 */
template<typename T>
class IMartyWilsonAdapter {
public:
    virtual ~IMartyWilsonAdapter() = default;

    /**
     * @brief Triggers a full calculation for a given Wilson basis and model.
     *
     * The exact pipeline (generation, compilation, numeric libs) is left to
     * concrete implementations.
     *
     * @param wilson      Wilson basis name.
     * @param model       Model name.
     * @param Q_match     Matching scale used in the numeric library.
     * @param model_path  Path to the model header used for template adaptation.
     */
    virtual void calculate(std::string wilson, std::string model, double Q_match, std::string model_path) = 0;

    /**
     * @brief Triggers a calculation where output label and model class differ.
     */
    virtual void calculate(std::string wilson,
                           std::string output_model,
                           std::string target_model,
                           double Q_match,
                           std::string model_path,
                           bool sm_like_filter,
                           bool bsm_split_generation = false,
                           bool full_target_generation = false) = 0;

    /**
     * @brief Returns the set of block names requiring special handling.
     *
     * @return Set of block names (same semantics as
     *         ::MartyInterface::get_special_blocks()).
     */
    virtual std::set<std::string>  get_special_blocks() = 0;

    /**
     * @brief Returns the set of dependencies for a given Wilson basis.
     *
     * These dependencies typically correspond to the parameters used in
     * the numeric MARTY example (e.g. ::InterpretedParam).
     *
     * @param wilson Wilson basis name.
     * @return Unordered set of dependencies of type @p T.
     */
    virtual std::unordered_set<T> get_dependencies(std::string wilson) = 0;


};

#endif