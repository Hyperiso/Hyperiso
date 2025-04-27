#ifndef IPARAMMODIFIER_H
#define IPARAMMODIFIER_H

/**
 * @example param_mutation_example.cpp
 * @brief Example usage of ParameterSetter and ParameterShifter
 * @defgroup DataMutationModule Data Mutation System
 * @brief Provides interfaces and classes to modify or shift parameter values dynamically.
 *
 * This module defines:
 * - Abstract interfaces for modifying parameters.
 * - Concrete implementations for setting or shifting parameters during runtime.
 *
 * ## Related Classes
 * - @ref IDataMutator
 * - @ref ParameterSetter
 * - @ref ParameterShifter
 * ## Diagram
 * @dot
 * digraph DataMutation {
 *   node [shape=record, fontname=Helvetica, fontsize=10];
 *
 *   IDataMutator [label="{ IDataMutator<T,U,V> | mutate(), change_mode() }"];
 *   ParameterSetter [label="{ ParameterSetter | mutate(), change_mode() }"];
 *   ParameterShifter [label="{ ParameterShifter | mutate(), change_mode() }"];
 *
 *   IDataMutator -> ParameterSetter;
 *   IDataMutator -> ParameterShifter;
 * }
 * @enddot
 */

 /**
 * @class IDataMutator
 * @ingroup DataMutationModule
 * @brief Interface for mutating parameters or changing their mode.
 *
 * @tparam T The parameter ID type.
 * @tparam U The type of the value to set (e.g., scalar).
 * @tparam V The type of the mode to apply (e.g., ParameterMode).
 */
template <typename T, typename U, typename V>
class IDataMutator {
public:
    virtual ~IDataMutator() = default;

    /**
     * @brief Mutates the value of a parameter.
     * @param param_id The identifier of the parameter to modify.
     * @param new_value The new value to set.
     */
    virtual void mutate(const T&, U) = 0;

    /**
     * @brief Changes the operational mode of a parameter.
     * @param param_id The identifier of the parameter.
     * @param new_mode The new mode to assign.
     */
    virtual void change_mode(const T&, V) = 0;
};

#endif // IPARAMMODIFIER_H
