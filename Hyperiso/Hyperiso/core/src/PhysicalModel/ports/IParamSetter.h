#ifndef __IPARAMSETTER_H__
#define __IPARAMSETTER_H__

#include "ParameterSetter.h"

class IParamSetter {
public:
    virtual ~IParamSetter() = default;

    virtual void set(double value) = 0;
};

#endif // __IPARAMSETTER_H__
