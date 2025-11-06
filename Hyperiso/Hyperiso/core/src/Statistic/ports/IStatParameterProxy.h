#ifndef ISTAT_PARAMETER_PROXY_H
#define ISTAT_PARAMETER_PROXY_H

#include "Math.h"
#include "ParameterProvider.h"

template<typename T, typename V>
class IStatParameterProxy {
public:
    template<typename... Args>
    requires HasCallableOperator<T, Args...>
    double operator()(Args&&... args) {
        return static_cast<T*>(this)->operator()(std::forward<Args>(args)...);
}
};

#endif