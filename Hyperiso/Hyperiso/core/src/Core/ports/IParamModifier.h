#ifndef __IPARAMMODIFIER_H__
#define __IPARAMMODIFIER_H__

#include "General.h"

class IParamModifier {
public:
    virtual ~IParamModifier() = default;

    virtual void modify(const ParamId&, double) = 0;
};

#endif // __IPARAMMODIFIER_H__
