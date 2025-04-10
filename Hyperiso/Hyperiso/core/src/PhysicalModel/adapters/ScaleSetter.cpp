#include "ScaleSetter.h"

void MatchingScaleSetter::set(double value) {
    ParameterSetter().mutate(ParamId{ParameterType::WILSON, "EW_SCALE", 1}, value);
}

void HadronicScaleSetter::set(double value) {
    ParameterSetter().mutate(ParamId{ParameterType::WILSON, "B_SCALE", 1}, value);
}