#include "Parameter.h"
#include <stdexcept>

Parameter::Parameter(ParamId id, double mean, double std_stat, double std_syst) 
    : id(id), expected(mean), deviation_stat(std_stat), deviation_syst(std_syst), value(mean), mode(ParameterMode::FIXED) {}


void Parameter::set_mode(ParameterMode new_mode) {
    if (new_mode != mode) {
        mode = new_mode;
        if (new_mode == ParameterMode::FIXED)
            value = expected;
    }
}
void Parameter::set_std(double stat, double syst) {
    this->deviation_stat = stat;
    this->deviation_syst = syst;
}


double Parameter::get_val() const {
    return value;
}

double Parameter::get_std() const {
    return std::hypot(deviation_stat, deviation_syst);
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
