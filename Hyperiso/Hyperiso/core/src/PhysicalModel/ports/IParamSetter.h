#ifndef __IPARAMSETTER_H__
#define __IPARAMSETTER_H__

#include "ParameterSetter.h"

template<typename T>
class IParamSetter {
public:
    virtual ~IParamSetter() = default;

    virtual void set(double value) = 0;
    virtual void switch_param(T scale_type) = 0;

protected:
    T param;
};

#endif // __IPARAMSETTER_H__
