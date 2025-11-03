#ifndef IPARAMSETTER_H
#define IPARAMSETTER_H

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

#endif // IPARAMSETTER_H
