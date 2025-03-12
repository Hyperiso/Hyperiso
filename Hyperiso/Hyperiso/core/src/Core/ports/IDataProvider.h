#ifndef __IPARAMETERPROVIDER_H__
#define __IPARAMETERPROVIDER_H__

#include <unordered_set>
#include <string>
#include <concepts>
#include <type_traits>
#include "General.h"
#include <utility>

template <typename, typename = std::void_t<>>
struct has_callable_operator : std::false_type {};

template <typename T>
struct has_callable_operator<T, std::void_t<decltype(std::declval<T>()(std::declval<int>()))>>
    : std::is_convertible<decltype(std::declval<T>()(std::declval<int>())), double> {};

template <typename T>
concept HasCallableOperator = has_callable_operator<T>::value;

template<typename T>
class IDataProvider {
public:
    virtual ~IDataProvider() = default;

    template<typename... Args>
    requires HasCallableOperator<T>
    double operator()(Args&&... args) {
        return static_cast<T*>(this)->operator()(std::forward<Args>(args)...);
    }
};


#endif // __IPARAMETERPROVIDER_H__
