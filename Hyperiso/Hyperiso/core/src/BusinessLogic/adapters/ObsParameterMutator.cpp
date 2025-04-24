#include "ObsParameterMutator.h"

void ObsParameterMutator::set(const ParamId &pid, scalar_t value) {
    this->p_set.mutate(pid, value);
}

void ObsParameterMutator::shift(const ParamId &pid, scalar_t value) {
    this->p_shift.mutate(pid, value);
}

void ObsParameterMutator::change_mode(const ParamId &pid, ParameterMode mod) {
    this->p_shift.change_mode(pid, mod);
}
