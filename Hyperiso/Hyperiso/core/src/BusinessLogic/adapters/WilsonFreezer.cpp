#include "WilsonFreezer.h"

void WilsonFreezer::freeze(WGroup group) {
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_EW_SCALE");
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_STANDARD");
    if (group == WGroup::B) {
        Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_TRADITIONAL");
    }    
}

void WilsonFreezer::unfreeze(WGroup group) {
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_EW_SCALE");
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_STANDARD");
    if (group == WGroup::B) {
        Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_B_SCALE_TRADITIONAL");
    }
}
