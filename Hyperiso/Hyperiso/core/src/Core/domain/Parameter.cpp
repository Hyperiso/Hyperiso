#include "Parameter.h"
#include "Block.h"
#include <stdexcept>
#include <utility>

Parameter::Parameter(ParamId id, scalar_t mean, scalar_t std_stat, scalar_t std_syst) 
    : id(id), expected(mean), deviation_stat(std_stat), deviation_syst(std_syst), shift(0), mode(ParameterMode::FIXED) {}


void Parameter::set_mode(ParameterMode new_mode) {
    if (new_mode != mode) {
        mode = new_mode;
    }
}

void Parameter::set_std(scalar_t stat, scalar_t syst) {
    this->deviation_stat = stat;
    this->deviation_syst = syst;
}

void Parameter::set_scale(double scale) {
    this->scale.emplace(scale);
}

void Parameter::set_bin(std::pair<double, double> bin) {
    this->binning.emplace(bin);
}

scalar_t Parameter::get_val() const {
    return this->mode == ParameterMode::FIXED ? expected : expected + shift;
}

void Parameter::set_expected(scalar_t val) {
    LOG_DEBUG("Parameter::set_expected of ", id.block, " ", id.code);
    this->expected = val;
    notifyObservers();
}

void Parameter::set_id(ParamId id) {
    this-> id = id;
}

scalar_t Parameter::get_combined_std() const {
    return std::hypot(static_cast<double>(deviation_stat),
                  static_cast<double>(deviation_syst));
}

std::pair<scalar_t, scalar_t> Parameter::get_std() const {
    return std::make_pair(deviation_stat, deviation_syst);
}

double Parameter::get_scale() {
    return this->scale.value_or(-1.0);
}

std::pair<double, double> Parameter::get_bin() {
    return this->binning.value_or(std::pair(-1.0, -1.0));
}

ParamId Parameter::get_id() const
{
    return id;
}

void Parameter::set_owner(ParameterType type) {
    id.set_parameter_type(type);
}

void Parameter::set_shift(scalar_t shift) {
    if (mode == ParameterMode::SHIFTABLE) {
        this->shift = shift;
    } else {
        throw std::runtime_error("Cannot change the value of immutable parameter " + id.block + " " + std::to_string(id.code));
    }
}

void Parameter::clear_below() {
    auto dependents = std::exchange(observers, {});
    LOG_DEBUG("cleaning param : ",this->id.block,this->id.code);

    this->clear_above();

    for (auto& obs : dependents) {
        if (obs) obs->clear_below();
    }
}



Parameter& Parameter::operator=(const Parameter& other) {
    this->id = other.id;
    this->set_expected(other.expected);
    this->deviation_stat = other.deviation_stat;
    this->deviation_syst = other.deviation_syst;
    this->mode = other.mode;
    this->shift = other.shift;
    this->scale = other.scale;
    this->binning = other.binning;
    return *this;
}

void Parameter::notifyObservers() {
    auto snapshot = observers;
    for (auto& observer : snapshot) {
        if (!observer) continue;
        LOG_DEBUG("Notifying observer", observer->id.block, observer->id.code,
                  "from parameter", id.block, id.code);
        observer->update();
    }
    observers.erase(std::remove(observers.begin(), observers.end(), nullptr), observers.end());
}


void Parameter::removeObserver(std::shared_ptr<Parameter> observer) { 
    auto it = std::find_if(observers.begin(), observers.end(),
        [&](const std::shared_ptr<Parameter>& p){ return p.get() == observer.get(); });
    if (it != observers.end()) observers.erase(it);
}


Parameter& Parameter::operator+=(const Parameter& other) {
    this->expected       += other.expected;
    this->deviation_stat = std::hypot(other.deviation_stat, this->deviation_stat);
    this->deviation_syst = std::hypot(other.deviation_syst, this->deviation_syst);
    this->shift          += other.shift;
    return *this;
}

Parameter &Parameter::operator*=(const scalar_t &scale) {
    this->expected       *= scale;
    this->deviation_stat *= std::abs(scale);
    this->deviation_syst *= std::abs(scale);
    this->shift          *= scale;
    return *this;
}
