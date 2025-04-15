#ifndef __IPARAMETERPROVIDER_H__
#define __IPARAMETERPROVIDER_H__

#include <unordered_set>
#include <string>
#include <concepts>
#include <type_traits>
#include "General.h"
#include "scalar.h"
#include <utility>

template <typename, typename... Args>
struct has_callable_operator : std::false_type {};

template <typename T, typename... Args>
struct has_callable_operator<T, std::void_t<decltype(std::declval<T>()(std::declval<Args>()...))>, Args...>
    : std::is_convertible<decltype(std::declval<T>()(std::declval<Args>()...)), scalar_t> {};

template <typename T, typename... Args>
concept HasCallableOperator = has_callable_operator<T,void, Args...>::value;

template <typename, typename... Args>
struct has_exists_function : std::false_type {};

template <typename T, typename... Args>
struct has_exists_function<T, std::void_t<decltype(std::declval<T>().exist(std::declval<Args>()...))>, Args...>
    : std::is_convertible<decltype(std::declval<T>().exist(std::declval<Args>()...)), bool> {};

template <typename T, typename... Args>
concept HasExistsFunction = has_exists_function<T, void, Args...>::value;

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


#endif // __IPARAMETERPROVIDER_H__
