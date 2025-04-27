#ifndef IPARAMETERPROVIDER_H
#define IPARAMETERPROVIDER_H

#include <unordered_set>
#include <string>
#include <concepts>
#include <type_traits>
#include "General.h"
#include "scalar.h"
#include <utility>

/**
 * @example data_provider_example.cpp
 * @brief Example usage of ParameterProvider, CorrelationProvider, and QCDProvider.
 * @defgroup DataProvidersModule Data Providers System
 * @brief Provides flexible access to parameter, correlation, and QCD data.
 *
 * This module defines a unified way to access different types of physical data
 * (parameters, correlations, QCD constants) through standardized callable interfaces.
 *
 * ## Overview
 *
 * - **IDataProvider** is a generic interface offering a call operator and optional existence checks.
 * - **ParameterProvider** retrieves parameter values from Parameters instances.
 * - **CorrelationProvider** retrieves correlations between parameters or observables.
 * - **QCDProvider** provides strong coupling constants and running masses.
 *
 * ## Related Classes
 * - @ref IDataProvider
 * - @ref ParameterProvider
 * - @ref CorrelationProvider
 * - @ref QCDProvider
 * - @ref IQCDProvider
 * ## Diagram
 * @dot
 * digraph DataProviders {
 *   node [shape=record, fontname=Helvetica, fontsize=10];
 *
 *   IDataProvider [label="{ IDataProvider<T> | ( ) , exists() }"];
 *   ParameterProvider [label="{ ParameterProvider | ( ) , exists() }"];
 *   CorrelationProvider [label="{ CorrelationProvider | ( ) }"];
 *   QCDProvider [label="{ QCDProvider | ( ) }"];
 *   IQCDProvider [label="{ IQCDProvider | get_constants() }"];
 *
 *   IDataProvider -> ParameterProvider;
 *   IDataProvider -> CorrelationProvider;
 *   IDataProvider -> QCDProvider;
 *   IQCDProvider -> QCDProvider;
 * }
 * @enddot
 */

template <typename, typename... Args>
struct has_callable_operator : std::false_type {};

/**
 * @brief Calls the underlying data provider with the given arguments.
 * @tparam Args Arguments to pass to the call operator.
 * @param args Arguments forwarded to the provider.
 * @return The value provided by the underlying provider.
 */
template <typename T, typename... Args>
struct has_callable_operator<T, std::void_t<decltype(std::declval<T>()(std::declval<Args>()...))>, Args...>
    : std::is_convertible<decltype(std::declval<T>()(std::declval<Args>()...)), scalar_t> {};

template <typename T, typename... Args>
concept HasCallableOperator = has_callable_operator<T,void, Args...>::value;

template <typename, typename... Args>
struct has_exists_function : std::false_type {};

/**
 * @brief Checks whether a given entry exists in the data provider.
 * @tparam Args Arguments to pass to the exists method.
 * @param args Arguments forwarded to the provider.
 * @return True if the entry exists, false otherwise.
 */
template <typename T, typename... Args>
struct has_exists_function<T, std::void_t<decltype(std::declval<T>().exist(std::declval<Args>()...))>, Args...>
    : std::is_convertible<decltype(std::declval<T>().exist(std::declval<Args>()...)), bool> {};

template <typename T, typename... Args>
concept HasExistsFunction = has_exists_function<T, void, Args...>::value;

/**
 * @class IDataProvider
 * @ingroup DataProvidersModule
 * @brief Generic interface for providing data through callable and exist methods.
 *
 * @tparam T Type implementing the actual data access logic.
 */
template<typename T>
class IDataProvider {
public:
    virtual ~IDataProvider() = default;

    template<typename... Args>
    requires HasCallableOperator<T, Args...>
    double operator()(Args&&... args) {
        return static_cast<T*>(this)->operator()(std::forward<Args>(args)...);
    }

    template<typename... Args>
    requires HasExistsFunction<T, Args...>
    bool exists(Args&&... args) const {
        return static_cast<T*>(this)->exists(std::forward<Args>(args)...);
    }
};


#endif // IPARAMETERPROVIDER_H
