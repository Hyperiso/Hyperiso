#include "Observable.h"

scalar_t Observable::eval() const {
    decay_parent->enable();
    return decay_parent->compute_observable(id);
}