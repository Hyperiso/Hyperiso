#include "WilsonFreezer.h"

void WilsonFreezer::freeze(WGroup group) {
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::MATCHING));
    for (WilsonBasis basis : w_proxy->get_bases(group)) {
        Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::HADRONIC, basis));
    }
}

void WilsonFreezer::unfreeze(WGroup group) {
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::MATCHING));
    for (WilsonBasis basis : w_proxy->get_bases(group)) {
        Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::HADRONIC, basis));
    }
}
