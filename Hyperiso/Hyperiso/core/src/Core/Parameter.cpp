#include "Parameter.h"
#include <stdexcept>

Parameter::Parameter(ParamId id, double mean, double std) 
    : id(id), expected(mean), deviation(std), value(mean), mode(ParameterMode::FIXED) {}


void Parameter::set_mode(ParameterMode new_mode) {
    if (new_mode != mode) {
        mode = new_mode;
        if (new_mode == ParameterMode::FIXED)
            value = expected;
    }
} 

double Parameter::get_val() const {
    return value;
}

double Parameter::get_std() const {
    return deviation;
}

ParamId Parameter::get_id() const {
    return id;
}

void Parameter::shift(double shift) {
    if (mode == ParameterMode::SHIFTABLE) {
        value = expected + shift;
    } else {
        throw std::runtime_error("Cannot change the value of immutable parameter " + id.block + " " + std::to_string(id.code));
    }
}
