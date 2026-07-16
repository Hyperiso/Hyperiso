#ifndef PARAM_OPTIMIZER_ADAPTER_H
#define PARAM_OPTIMIZER_ADAPTER_H

#include "ParamOptimizer.h"
#include "IParamOptimizer.h"

/**
 * @file ParamOptimizerAdapter.h
 * @brief Adapter between IParamOptimizer and ParamOptimizer.
 *
 * This header declares ParamOptimizerAdapter, a concrete implementation of
 * IParamOptimizer that internally wraps a ParamOptimizer instance built from
 * a list of ParameterType scopes.
 */

 /**
 * @class ParamOptimizerAdapter
 * @ingroup ParamOptimizationModule
 * @brief Adapter that builds a ParamOptimizer from Parameters instances.
 *
 * ParamOptimizerAdapter:
 *  - takes a list of ParameterType scopes in its constructor,
 *  - retrieves the corresponding BlockAccessor instances from Parameters,
 *  - constructs a ParamOptimizer working over those scopes,
 *  - exposes the IParamOptimizer interface for staged parameter updates.
 *
 * Example:
 * @code
 * ParamOptimizerAdapter optimizer({ParameterType::SM, ParameterType::BSM});
 * optimizer.set_value("MASS", LhaID(25), 123.4);
 * optimizer.commit(); // propagate changes in a batched way
 * @endcode
 */
class ParamOptimizerAdapter : public IParamOptimizer {
public:
    /**
     * @brief Constructs the adapter from a list of parameter type scopes.
     *
     * For each ParameterType in @p scopes, the constructor retrieves
     * Parameters::GetInstance(scope)->get_block_accessor() and passes the
     * resulting collection of BlockAccessor pointers to an internal
     * ParamOptimizer instance.
     *
     * @param scopes List of ParameterTypes to include in the optimization scope.
     */
    ParamOptimizerAdapter(std::vector<ParameterType> scopes);
    
    /// @copydoc IParamOptimizer::set_value()
    void set_value(const BlockName& block, const LhaID& id, scalar_t v) override;

    /// @copydoc IParamOptimizer::set_param()
    void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) override;

    /// @copydoc IParamOptimizer::remove()
    void remove(const BlockName& block, const LhaID& id) override;

    /// @copydoc IParamOptimizer::commit()
    void commit(bool coalesce = true) override;

    /// @copydoc IParamOptimizer::clear()
    void clear() override;

private:
    /// Underlying concrete optimizer operating on BlockAccessor scopes.
    std::shared_ptr<ParamOptimizer> po;

};

#endif