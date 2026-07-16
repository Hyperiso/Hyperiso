#include "WilsonFreezer.h"

void WilsonFreezer::freeze(WGroupId group) {
    Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::MATCHING));
    for (WilsonBasis basis : w_proxy->get_bases(group)) {
        Freezer::freeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::HADRONIC, basis));
    }
}

void WilsonFreezer::unfreeze(WGroupId group) {
    Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::MATCHING));
    for (WilsonBasis basis : w_proxy->get_bases(group)) {
        Freezer::unfreeze(ParameterType::WILSON, GroupMapper::str(group, ScaleType::HADRONIC, basis));
    }
}
