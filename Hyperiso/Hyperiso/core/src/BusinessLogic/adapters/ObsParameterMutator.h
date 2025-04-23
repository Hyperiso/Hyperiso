#include "IObsParameterMutator.h"
#include "ParameterSetter.h"


class ObsParameterMutator: public IObsParameterMutator {
public:
void mutate(const ParamId&, double) override;
};