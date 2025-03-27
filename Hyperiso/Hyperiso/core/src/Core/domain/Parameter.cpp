#include "Parameter.h"
#include <stdexcept>

Parameter::Parameter(ParamId id, double mean, double std_stat, double std_syst) 
    : id(id), expected(mean), deviation_stat(std_stat), deviation_syst(std_syst), shift(0), mode(ParameterMode::FIXED) {}


void Parameter::set_mode(ParameterMode new_mode) {
    if (new_mode != mode) {
        mode = new_mode;
    }
}

void Parameter::set_std(double stat, double syst) {
    this->deviation_stat = stat;
    this->deviation_syst = syst;
}

double Parameter::get_val() const {
    return this->mode == ParameterMode::FIXED ? expected : expected + shift;
}

void Parameter::set_expected(double val) {
    this->expected = val;
    notifyObservers();
}

double Parameter::get_std() const {
    return std::hypot(deviation_stat, deviation_syst);
}

ParamId Parameter::get_id() const {
    return id;
}

void Parameter::set_owner(ParameterType type) {
    id.set_parameter_type(type);
}

void Parameter::set_shift(double shift) {
    if (mode == ParameterMode::SHIFTABLE) {
        this->shift = shift;
    } else {
        throw std::runtime_error("Cannot change the value of immutable parameter " + id.block + " " + std::to_string(id.code));
    }
}
