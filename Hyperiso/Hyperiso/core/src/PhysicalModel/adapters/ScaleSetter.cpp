#include "ScaleSetter.h"

void ScaleSetter::set(double value) {
    ParamId pid {ParameterType::WILSON, ScaleTypeMapper::block(this->param), 1};
    ParameterSetter().mutate(pid, value);
}