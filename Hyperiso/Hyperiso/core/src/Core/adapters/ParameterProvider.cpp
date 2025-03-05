#include "ParameterProvider.h"

double ParameterProvider::operator()(const ParamId &pid, DataType type) {
    if (this->p_type) {
        LOG_WARN("LogicError", "This ParameterProvider already has a type.");
    }

}

double ParameterProvider::operator()(const std::string &block, const LhaID &id, DataType type) {
    if (!this->p_type) {
        LOG_ERROR("LogicError", "Please specify a parameter type for the ParameterProvider.");
    }

    return 0.0;
}
