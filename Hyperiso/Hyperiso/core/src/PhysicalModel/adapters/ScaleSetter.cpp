#include "ScaleSetter.h"

void ScaleSetter::set(double value) {
    ParameterSetter().mutate(ParamId{ParameterType::WILSON, ScaleBlockMapper::block(this->scale_type), 1}, value);
}