#include "WilsonFreezer.h"

//TODO or NOT TODO : need to had traditionnal. but do we need full as well ?
void WilsonFreezer::freeze(WGroup group) {
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_EW_SCALE");
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_STANDARD");
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_TRADITIONAL");
}

void WilsonFreezer::unfreeze(WGroup group) {
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_EW_SCALE");
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_STANDARD");
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_TRADITIONAL");
}
