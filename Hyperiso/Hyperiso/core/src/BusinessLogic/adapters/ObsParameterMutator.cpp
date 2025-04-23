#include "ObsParameterMutator.h"

void ObsParameterMutator::mutate(const ParamId& pid, double value) {
    ParameterSetter().mutate(pid, value);
}