#ifndef IPARAM_ADAPTER_H
#define IPARAM_ADAPTER_H

template<typename T, typename V>
class IParamAdapter {
public:
    virtual double operator()(T x, V y) = 0;
};

#endif