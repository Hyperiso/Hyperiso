#ifndef IDATABASe_ADAPTER_H
#define IDATABASe_ADAPTER_H

#include "scalar.h"
#include "ParameterProvider.h"

template<typename T, typename V>
class IDataBaseProxy {
public:
    template<typename... Args>
    requires HasCallableOperator<T, Args...>
    double operator()(Args&&... args) {
        return static_cast<T*>(this)->operator()(std::forward<Args>(args)...);
}
};

#endif