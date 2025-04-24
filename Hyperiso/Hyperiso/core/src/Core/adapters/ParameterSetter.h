#ifndef PARAMETERMODIFIER_H
#define PARAMETERMODIFIER_H

#include "IDataMutator.h"
#include "Parameters.h"
#include "Math.h"

class ParameterSetter : public IDataMutator<ParamId, scalar_t, ParameterMode> {
public:
    void mutate(const ParamId& pid, scalar_t value) override;
    void change_mode(const ParamId& pid, ParameterMode mode) override;
};

#endif // __PARAMETERMODIFIER_H__
