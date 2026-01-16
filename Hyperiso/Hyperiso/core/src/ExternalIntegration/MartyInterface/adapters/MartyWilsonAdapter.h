#ifndef MARTY_WILSON_ADAPTER_H
#define MARTY_WILSON_ADAPTER_H

#include "MartyInterface.h"
#include "IMartyWilsonAdapter.h"

/**
 * @file MartyWilsonAdapter.h
 * @brief Declares a concrete IMartyWilsonAdapter implementation built on MartyInterface.
 *
 * This adapter wraps ::MartyInterface and exposes it via the generic
 * ::IMartyWilsonAdapter interface, making it easy to swap implementations
 * or mock the adapter in tests.
 */

/**
 * @class MartyWilsonAdapter
 * @ingroup CodeGenerationModule
 * @brief Concrete adapter for Wilson calculations based on ::MartyInterface.
 *
 * This class simply delegates:
 *  - ::calculate to ::MartyInterface::calculate,
 *  - ::get_special_blocks to ::MartyInterface::get_special_blocks,
 *  - ::get_dependencies to ::MartyInterface::get_dependencies.
 *
 * It is parameterized with ::InterpretedParam in the ::IMartyWilsonAdapter
 * base template.
 */
class MartyWilsonAdapter : public IMartyWilsonAdapter<InterpretedParam> {
public:
    /**
     * @brief Default constructor creating an internal ::MartyInterface.
     *
     * Uses the default constructor of ::MartyInterface, which wires all
     * default dependencies (ModelAPI, proxies, ports).
     */
    MartyWilsonAdapter() {martyInterface = MartyInterface();}

    /**
     * @brief Launches a full calculation for a given (Wilson, model) pair.
     *
     * Delegates directly to ::MartyInterface::calculate.
     */
    void calculate(std::string wilson, std::string model, double Q_match, std::string model_path) override {martyInterface.calculate(wilson, model, Q_match, model_path);}

    /**
     * @brief Returns the set of special blocks requiring custom handling.
     *
     * Delegates to ::MartyInterface::get_special_blocks.
     */
    std::set<std::string>  get_special_blocks() override {return martyInterface.get_special_blocks();}

    /**
     * @brief Returns the dependencies associated with a given Wilson basis.
     *
     * Delegates to ::MartyInterface::get_dependencies.
     */
    std::unordered_set<InterpretedParam> get_dependencies(std::string wilson) override {return martyInterface.get_dependencies(wilson);}

private:
    /// Underlying façade used to perform all MARTY-related work.
    MartyInterface martyInterface;
};

#endif  // MARTY_WILSON_ADAPTER_H