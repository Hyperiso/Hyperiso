#ifndef __OBSPARAMETERMUTATOR_H__
#define __OBSPARAMETERMUTATOR_H__

#include "Include.h"
#include "Math.h"
#include "IObsParameterMutator.h"
#include "ParameterSetter.h"
#include "ParameterShifter.h"

class ObsParameterMutator: public IObsParameterMutator<ParamId, scalar_t> {
public:
    void set(const ParamId&, scalar_t) override;
    void shift(const ParamId&, scalar_t) override;

private:
    ParameterSetter p_set;
    ParameterShifter p_shift;
};

#endif // __OBSPARAMETERMUTATOR_H__