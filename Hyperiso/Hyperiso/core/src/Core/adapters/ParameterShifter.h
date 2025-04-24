#ifndef __PARAMETERSHIFTER_H__
#define __PARAMETERSHIFTER_H__

#include "IDataMutator.h"
#include "Parameters.h"

class ParameterShifter : public IDataMutator {
public:
    void mutate(const ParamId& pid, scalar_t value) override;
};


#endif // __PARAMETERSHIFTER_H__
