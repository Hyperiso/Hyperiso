#include "WilsonFreezer.h"

void WilsonFreezer::freeze(WGroup group) {
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_MATCH");
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group) + "_HADRONIC");
}

void WilsonFreezer::unfreeze(WGroup group) {
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_MATCH");
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group) + "_HADRONIC");
}
