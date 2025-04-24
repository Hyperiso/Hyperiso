#ifndef PARAMETERSHIFTER_H
#define PARAMETERSHIFTER_H

#include "IDataMutator.h"
#include "Parameters.h"

class ParameterShifter : public IDataMutator<ParamId, scalar_t, ParameterMode> {
public:
    void mutate(const ParamId& pid, scalar_t value) override;
    void change_mode(const ParamId& pid, ParameterMode mode) override;
};


#endif // __PARAMETERSHIFTER_H__
