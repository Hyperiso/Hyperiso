#ifndef MARTY_WILSON_PROXY_H
#define MARTY_WILSON_PROXY_H

#include "IMartyWilsonProxy.h"
#include "MartyWilsonAdapter.h"

#include <utility>

/**
 * @file MartyWilsonProxy.h
 * @brief Concrete proxy for Wilson coefficient calculations via MARTY.
 *
 * This header defines @ref MartyWilsonProxy, a thin façade implementing
 * @ref IMartyWilsonProxy and delegating all logic to
 * @ref MartyWilsonAdapter.
 */

/**
 * @class MartyWilsonProxy
 * @ingroup MartyIntegrationModule
 * @brief Concrete proxy forwarding Wilson calculations to MARTY.
 *
 * This class acts as a stable interface layer between higher-level
 * Hyperiso components and the MARTY-specific implementation.
 *
 * Responsibilities:
 *  - delegate calculation requests,
 *  - expose special blocks,
 *  - expose parameter dependencies.
 */
class MartyWilsonProxy : public IMartyWilsonProxy<InterpretedParam> {
public:
    /**
     * @brief Default constructor.
     *
     * Initializes the internal @ref MartyWilsonAdapter.
     */
    MartyWilsonProxy() {martyAdapter = MartyWilsonAdapter();}

    /// @copydoc IMartyWilsonProxy::calculate
    void calculate(std::string wilson, std::string model, double Q_match, std::string model_path) override {martyAdapter.calculate(wilson, model, Q_match, model_path);}

    /// @copydoc IMartyWilsonProxy::calculate(std::string,std::string,std::string,double,std::string,bool)
    void calculate(std::string wilson,
                   std::string output_model,
                   std::string target_model,
                   double Q_match,
                   std::string model_path,
                   bool sm_like_filter,
                   bool bsm_split_generation = false,
                   bool full_target_generation = false) override {
        martyAdapter.calculate(
            wilson,
            output_model,
            target_model,
            Q_match,
            model_path,
            sm_like_filter,
            bsm_split_generation,
            full_target_generation
        );
    }

    std::string calculate_isolated(std::string wilson,
                                   std::string output_model,
                                   std::string target_model,
                                   double Q_match,
                                   std::string model_path,
                                   bool sm_like_filter,
                                   bool bsm_split_generation = false,
                                   bool full_target_generation = false) override {
        return martyAdapter.calculate_isolated(
            std::move(wilson),
            std::move(output_model),
            std::move(target_model),
            Q_match,
            std::move(model_path),
            sm_like_filter,
            bsm_split_generation,
            full_target_generation
        );
    }

    /// @copydoc IMartyWilsonProxy::get_special_blocks
    std::set<std::string>  get_special_blocks() override {return martyAdapter.get_special_blocks();}

    /// @copydoc IMartyWilsonProxy::get_dependencies
    std::unordered_set<InterpretedParam> get_dependencies(std::string wilson) override {return martyAdapter.get_dependencies(wilson);}

private:
    /// Internal adapter implementing the actual MARTY logic.
    MartyWilsonAdapter martyAdapter;
};

#endif // MARTY_WILSON_PROXY_H