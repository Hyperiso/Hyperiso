#ifndef IPARAMMODIFIER_H
#define IPARAMMODIFIER_H

/**
 * @file IDataMutator.h
 * @brief Interface for mutating parameter values and modes.
 *
 * This header declares the templated interface ::IDataMutator, which defines
 * a generic protocol for:
 *  - changing the numerical value of a parameter-like entity,
 *  - changing its operational mode (for instance, fixed vs shiftable).
 *
 * It is meant to be specialized by concrete mutators such as:
 *  - ParameterSetter (direct assignment),
 *  - ParameterShifter (incremental shift).
 */

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
 * This is a generic interface that decouples the notion of “what” is mutated
 * from “how” the mutation is implemented.
 *
 * Typical template arguments are:
 *  - `T` : an identifier type (e.g. ParamId),
 *  - `U` : a value type (e.g. scalar_t),
 *  - `V` : a mode type (e.g. ParameterMode).
 *
 * Concrete implementations can then decide whether they:
 *  - set absolute values (ParameterSetter),
 *  - apply relative shifts (ParameterShifter),
 *  - or do more advanced logic (e.g. bounds checking, logging).
 *
 * @tparam T The parameter ID type.
 * @tparam U The type of the value to set (e.g., scalar).
 * @tparam V The type of the mode to apply (e.g., ParameterMode).
 */
template <typename T, typename U, typename V>
class IDataMutator {
public:
    /// Virtual destructor for interface.
    virtual ~IDataMutator() = default;

    /**
     * @brief Mutates the value of a parameter.
     *
     * Implementations are free to interpret @p new_value as:
     *  - an absolute replacement (e.g. set to new_value),
     *  - a relative shift (e.g. add new_value),
     *  - or something more elaborate, depending on their semantics.
     *
     * @param param_id  The identifier of the parameter to modify.
     * @param new_value The new value to apply (interpretation is implementation-dependent).
     */
    virtual void mutate(const T&, U) = 0;

    /**
     * @brief Changes the operational mode of a parameter.
     *
     * Typical modes are represented by an enum (e.g. ParameterMode),
     * and control how the parameter behaves in fits, scans, or updates.
     *
     * @param param_id The identifier of the parameter.
     * @param new_mode The new mode to assign.
     */
    virtual void change_mode(const T&, V) = 0;
};

#endif // IPARAMMODIFIER_H
