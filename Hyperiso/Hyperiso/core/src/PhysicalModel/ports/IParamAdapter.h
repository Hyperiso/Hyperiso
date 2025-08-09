#ifndef IPARAM_ADAPTER_H
#define IPARAM_ADAPTER_H
#include "scalar.h"

template<typename T, typename V>
class IParameterProxy {
public:
    virtual scalar_t operator()(const T& x, const V& y) const = 0;
    virtual bool exist(const T& block, const V& id) const = 0;
    virtual double get_scale(const T& block) const = 0;
};

#endif