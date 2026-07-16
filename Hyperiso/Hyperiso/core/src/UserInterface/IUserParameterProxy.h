#ifndef IUSER_PARAM_PROXY_H
#define IUSER_PARAM_PROXY_H

#include "scalar.h"

template<typename T, typename V>
class IUserParameterProxy {
public:
    virtual std::optional<double> get_value(T block, V id) = 0;

    virtual void set_value(T block, V id, double val) = 0;
};

#endif // IPARAM_PROXY_H