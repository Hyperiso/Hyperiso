#ifndef __PARAMETERMODIFIER_H__
#define __PARAMETERMODIFIER_H__

#include "IDataMutator.h"
#include "Parameters.h"

class ParameterSetter : public IDataMutator {
public:
    void mutate(const ParamId& pid, double value);
};

#endif // __PARAMETERMODIFIER_H__
