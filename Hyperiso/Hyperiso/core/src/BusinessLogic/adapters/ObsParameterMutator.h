#ifndef OBSPARAMETERMUTATOR_H
#define OBSPARAMETERMUTATOR_H

#include "Include.h"
#include "Math.h"
#include "IObsParameterMutator.h"
#include "ParameterSetter.h"
#include "ParameterShifter.h"

class ObsParameterMutator: public IObsParameterMutator<ParamId, scalar_t, ParameterMode> {
public:
    void set(const ParamId&, scalar_t) override;
    void shift(const ParamId&, scalar_t) override;
    void change_mode(const ParamId&, ParameterMode) override;

private:
    ParameterSetter p_set;
    ParameterShifter p_shift;
};

#endif // __OBSPARAMETERMUTATOR_H__