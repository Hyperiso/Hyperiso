#ifndef __IPARAMETERPROVIDER_H__
#define __IPARAMETERPROVIDER_H__

#include <unordered_set>
#include <string>
#include <concepts>
#include <type_traits>
#include "General.h"

template<typename T, typename... Args>
concept HasCallableOperator = requires(T t, Args... args) {
    { t(args...) } -> std::convertible_to<double>;
};

template<typename T>
requires HasCallableOperator<T>
class IDataProvider {
public:
    virtual ~IDataProvider() = default;

    template<typename... Args>
    double operator()(Args... args) {
        return static_cast<T*>(this)->operator()(args...);
    }
};


// virtual std::unordered_set<std::string> get_block_list() = 0;
// virtual std::map<LhaID, double> get_block_values(const std::string&) = 0;
// virtual std::map<LhaID, Parameter> get_block_parameters(const std::string&) = 0;


#endif // __IPARAMETERPROVIDER_H__
