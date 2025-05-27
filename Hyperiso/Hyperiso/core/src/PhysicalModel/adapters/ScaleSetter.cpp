#include "ScaleSetter.h"

void ScaleSetter::set(double value) {
    LOG_INFO("ScaleSetter::set");
    ParamId pid {ParameterType::WILSON, ScaleTypeMapper::block(this->scale_type), 1};
    ParameterSetter().mutate(pid, value);
}