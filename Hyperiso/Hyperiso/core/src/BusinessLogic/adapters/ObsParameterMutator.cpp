#include "ObsParameterMutator.h"

void ObsParameterMutator::set(const ParamId &pid, scalar_t value) {
    this->p_set.mutate(pid, value);
}

void ObsParameterMutator::shift(const ParamId &pid, scalar_t value) {
    this->p_shift.mutate(pid, value);
}
