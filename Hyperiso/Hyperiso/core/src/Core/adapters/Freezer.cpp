#include "Freezer.h"

void Freezer::freeze(const ParameterType &p_type, const std::string &block_name) {
    Parameters::GetInstance(p_type)->freeze_block(block_name);
}

void Freezer::freeze(const ParamId &pid) {
    Parameters::GetInstance(pid.type.value())->freeze_param(pid.block, pid.code);
}

void Freezer::unfreeze(const ParameterType &p_type, const std::string &block_name) {
    Parameters::GetInstance(p_type)->unfreeze_block(block_name);
}

void Freezer::unfreeze(const ParamId &pid) {
    Parameters::GetInstance(pid.type.value())->unfreeze_param(pid.block, pid.code);
}
