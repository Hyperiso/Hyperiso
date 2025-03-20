#ifndef IPARAM_ADAPTER_H
#define IPARAM_ADAPTER_H

template<typename T, typename V>
class IParameterProxy {
public:
    virtual double operator()(const T& x, const V& y) = 0;
};

#endif