#ifndef __IPARAMMODIFIER_H__
#define __IPARAMMODIFIER_H__

#include "General.h"

class IDataMutator {
public:
    virtual ~IDataMutator() = default;

    virtual void mutate(const ParamId&, scalar_t) = 0;
};

#endif // __IPARAMMODIFIER_H__
