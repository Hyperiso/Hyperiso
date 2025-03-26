#ifndef IPARAM_ADAPTER_H
#define IPARAM_ADAPTER_H

template<typename T, typename V>
class IParameterProxy {
public:
    virtual double operator()(const T& x, const V& y) const = 0 ;
};

#endif