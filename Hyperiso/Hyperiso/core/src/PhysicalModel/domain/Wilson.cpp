#include "Wilson.h"

void WilsonCoefficient::set_owned(bool owned) {
    if (this->is_owned && owned) {
        LOG_ERROR("LogicError", "WilsonCoefficient is already owned by a WilsonGroup and cannot be shared.");
    }

    this->is_owned = owned;
}
