#ifndef IOBS_PARAMETER_MUTATOR_H
#define IOBS_PARAMETER_MUTATOR_H

#include "General.h"

class IObsParameterMutator {
public:
    virtual ~IObsParameterMutator() = default;

    virtual void mutate(const ParamId&, double) = 0;
};

#endif // __IPARAMMODIFIER_H__
